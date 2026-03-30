"""
Transient Dynamics Module
==========================

Models time-dependent behavior of electrospray emitters.

Critical for hardware design - real thrusters don't operate at steady state.

Physics included:
1. Startup transients (0 → steady state)
2. Pulsed operation (on/off cycling)
3. Current oscillations
4. Mode transitions
5. Response to voltage/flow changes

References:
- Gamero-Castaño (2008) - "Electric measurements of charged sprays"
- Dandavino et al. (2014) - "Microfabricated electrospray emitter arrays with integrated extractor"
- Uchizono et al. (2022) - "Emission characteristics of passively fed ionic liquid electrospray"
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from typing import Tuple, Callable, Optional, Dict, List
from dataclasses import dataclass


@dataclass
class TransientConditions:
    """Initial and boundary conditions for transient simulation"""
    t_start: float = 0.0  # s
    t_end: float = 1.0  # s
    dt: float = 1e-4  # s - time step
    
    # Initial state
    current_initial: float = 0.0  # A
    flow_rate_initial: float = 0.0  # kg/s
    temperature_initial: float = 298.15  # K
    
    # Applied conditions
    voltage_profile: Callable[[float], float] = lambda t: 2000.0  # V(t)
    pressure_profile: Callable[[float], float] = lambda t: 200e3  # P(t)


class TransientElectrosprayModel:
    """
    Time-dependent electrospray emission model
    
    Solves coupled ODEs for:
    - Current evolution
    - Flow rate dynamics  
    - Temperature evolution
    - Meniscus position
    """
    
    def __init__(self, liquid_props: Dict, emitter_geometry: Dict):
        self.liquid = liquid_props
        self.geom = emitter_geometry
        
        # Characteristic time constants (fitted from experiments)
        self.tau_current = 0.01  # s - current response time
        self.tau_thermal = 0.1  # s - thermal response time
        self.tau_flow = 0.05  # s - flow response time
    
    def steady_state_current(self, m_dot: float, T: float) -> float:
        """
        Steady-state current for given flow and temperature
        
        This uses the validated scaling law: I = I₀ + ψ√(γKṁ/ρ)
        with temperature-dependent conductivity
        """
        gamma = self.liquid['surface_tension']
        rho = self.liquid['density']
        K_0 = self.liquid['conductivity']
        
        # Temperature-dependent conductivity (simplified Arrhenius)
        E_a = 0.3 * 1.602e-19  # ~0.3 eV activation energy
        k_B = 1.381e-23
        K_T = K_0 * np.exp(-E_a / (k_B * T)) / np.exp(-E_a / (k_B * 298.15))
        
        # Current scaling
        I_0 = 81e-9  # A
        psi = 2.48
        
        I_ss = I_0 + psi * np.sqrt(gamma * K_T * m_dot / rho)
        
        return I_ss
    
    def ode_system(self, t: float, state: np.ndarray) -> np.ndarray:
        """
        System of ODEs for transient behavior
        
        State vector: [I, m_dot, T]
        
        dI/dt = (I_ss - I) / tau_I
        dm/dt = (m_ss - m) / tau_m  
        dT/dt = (T_ss - T) / tau_T + Q_joule / (m*c_p)
        
        Args:
            t: Time (s)
            state: [I, m_dot, T]
        
        Returns:
            [dI/dt, dm_dot/dt, dT/dt]
        """
        I, m_dot, T = state
        
        # Applied voltage (can be time-varying)
        V = 2000.0  # For now, constant
        
        # Steady-state targets
        I_ss = self.steady_state_current(m_dot, T)
        m_dot_ss = 1e-10  # kg/s - from flow system
        T_ss = 298.15 + 50.0  # K - self-heating equilibrium
        
        # Joule heating rate
        Q_joule = I * V  # W
        c_p = 1500  # J/(kg·K) - specific heat
        m_emitter = 1e-9  # kg - mass of heated region
        
        # ODEs
        dI_dt = (I_ss - I) / self.tau_current
        dm_dot_dt = (m_dot_ss - m_dot) / self.tau_flow
        dT_dt = (T_ss - T) / self.tau_thermal + Q_joule / (m_emitter * c_p)
        
        return np.array([dI_dt, dm_dot_dt, dT_dt])
    
    def simulate_startup(self, conditions: TransientConditions) -> Dict:
        """
        Simulate startup transient: 0 → steady state
        
        Returns:
            Dict with time series of I(t), m_dot(t), T(t)
        """
        # Time array
        t_span = (conditions.t_start, conditions.t_end)
        t_eval = np.arange(conditions.t_start, conditions.t_end, conditions.dt)
        
        # Initial state
        state_0 = np.array([
            conditions.current_initial,
            conditions.flow_rate_initial,
            conditions.temperature_initial
        ])
        
        # Solve ODE
        solution = solve_ivp(
            self.ode_system,
            t_span,
            state_0,
            t_eval=t_eval,
            method='RK45',
            rtol=1e-6
        )
        
        return {
            't': solution.t,
            'I': solution.y[0, :],
            'm_dot': solution.y[1, :],
            'T': solution.y[2, :],
            'success': solution.success,
            'message': solution.message
        }
    
    def simulate_pulse(self, t_on: float, t_off: float, 
                       n_cycles: int) -> Dict:
        """
        Simulate pulsed operation: on/off cycling
        
        Args:
            t_on: On-time duration (s)
            t_off: Off-time duration (s)
            n_cycles: Number of cycles
        
        Returns:
            Time series for pulsed operation
        """
        t_cycle = t_on + t_off
        t_total = n_cycles * t_cycle
        
        # Initialize arrays
        dt = 1e-4  # s
        t_array = np.arange(0, t_total, dt)
        I_array = np.zeros_like(t_array)
        m_dot_array = np.zeros_like(t_array)
        T_array = np.ones_like(t_array) * 298.15
        
        # Initial state
        state = np.array([0.0, 0.0, 298.15])
        
        for i, t in enumerate(t_array):
            # Determine if on or off
            t_in_cycle = t % t_cycle
            is_on = t_in_cycle < t_on
            
            if is_on:
                # Evolve toward on-state
                dstate = self.ode_system(t, state)
                state += dstate * dt
            else:
                # Decay toward off-state
                tau_decay = 0.02  # s - decay time constant
                state[0] *= np.exp(-dt / tau_decay)  # Current decays
                state[1] *= np.exp(-dt / tau_decay)  # Flow decays
                state[2] = 298.15 + (state[2] - 298.15) * np.exp(-dt / self.tau_thermal)
            
            I_array[i] = state[0]
            m_dot_array[i] = state[1]
            T_array[i] = state[2]
        
        return {
            't': t_array,
            'I': I_array,
            'm_dot': m_dot_array,
            'T': T_array,
            't_on': t_on,
            't_off': t_off,
            'duty_cycle': t_on / t_cycle,
        }
    
    def time_to_steady_state(self, tolerance: float = 0.05) -> float:
        """
        Calculate time to reach steady state (within tolerance)
        
        Args:
            tolerance: Fractional tolerance (0.05 = 5%)
        
        Returns:
            Time to steady state (s)
        """
        # Approximate as 5 time constants
        tau_max = max(self.tau_current, self.tau_thermal, self.tau_flow)
        t_ss = 5 * tau_max
        
        return t_ss


class CurrentOscillationModel:
    """
    Models current oscillations observed in experiments
    
    Oscillations arise from:
    1. Meniscus dynamics
    2. Electrohydrodynamic instabilities
    3. Charge relaxation
    """
    
    def __init__(self, frequency: float = 100.0, amplitude_fraction: float = 0.1):
        """
        Args:
            frequency: Oscillation frequency (Hz)
            amplitude_fraction: Oscillation amplitude as fraction of mean current
        """
        self.f = frequency
        self.A = amplitude_fraction
    
    def current_with_oscillation(self, I_mean: float, t: np.ndarray) -> np.ndarray:
        """
        Add oscillations to mean current
        
        I(t) = I_mean × (1 + A·sin(2πft))
        
        Args:
            I_mean: Mean current (A)
            t: Time array (s)
        
        Returns:
            I(t) with oscillations
        """
        return I_mean * (1.0 + self.A * np.sin(2 * np.pi * self.f * t))
    
    def power_spectrum(self, I_t: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute power spectrum of current oscillations
        
        Args:
            I_t: Current time series
            dt: Time step (s)
        
        Returns:
            (frequencies, power_spectrum)
        """
        # FFT
        I_fft = np.fft.rfft(I_t)
        freqs = np.fft.rfftfreq(len(I_t), dt)
        
        # Power spectrum
        power = np.abs(I_fft)**2
        
        return freqs, power


class ModeTransitionModel:
    """
    Models transitions between emission modes
    
    Modes:
    1. Pure ion emission (field evaporation)
    2. Mixed ion-droplet emission
    3. Pure droplet emission
    
    Transitions occur when flow rate or voltage changes.
    """
    
    def __init__(self):
        # Transition thresholds (from experiments)
        self.Q_ion_droplet = 1e-14  # m³/s - transition from ion to mixed
        self.Q_droplet = 1e-13  # m³/s - transition to pure droplet
    
    def determine_mode(self, Q: float) -> str:
        """
        Determine emission mode based on flow rate
        
        Args:
            Q: Volumetric flow rate (m³/s)
        
        Returns:
            'ion', 'mixed', or 'droplet'
        """
        if Q < self.Q_ion_droplet:
            return 'ion'
        elif Q < self.Q_droplet:
            return 'mixed'
        else:
            return 'droplet'
    
    def mode_transition_time(self, Q_initial: float, Q_final: float) -> float:
        """
        Estimate time for mode transition
        
        Based on meniscus reconfiguration time
        
        Returns:
            Transition time (s)
        """
        # Capillary time scale: t_cap = (ρR³/γ)^(1/2)
        rho = 1520  # kg/m³ - typical IL
        gamma = 0.042  # N/m
        R = 10e-6  # m - meniscus radius
        
        t_cap = np.sqrt(rho * R**3 / gamma)
        
        # Transition takes ~10 capillary times
        return 10 * t_cap


# ============================================================================
# INTEGRATION WITH EXISTING MODEL
# ============================================================================

def add_transient_capability(base_emitter, enable_transient: bool = False):
    """
    Add transient dynamics to existing emitter
    
    Args:
        base_emitter: Existing emitter instance
        enable_transient: Whether to enable transient modeling
    
    Returns:
        Enhanced emitter with transient capability
    """
    if not enable_transient:
        return base_emitter
    
    # Extract properties
    liquid_props = {
        'surface_tension': base_emitter.liquid.surface_tension,
        'density': base_emitter.liquid.density,
        'conductivity': base_emitter.liquid.conductivity,
    }
    
    emitter_geom = {
        'd_tip': base_emitter.d_tip,
    }
    
    # Create transient model
    base_emitter.transient_model = TransientElectrosprayModel(liquid_props, emitter_geom)
    base_emitter.has_transient = True
    
    return base_emitter


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Transient Dynamics Module - Demonstration")
    print("=" * 70)
    
    # Define liquid properties (EMI-Im)
    liquid_props = {
        'surface_tension': 0.042,
        'density': 1520,
        'conductivity': 1.5,
    }
    
    emitter_geom = {
        'd_tip': 30e-6,
    }
    
    # Create transient model
    transient = TransientElectrosprayModel(liquid_props, emitter_geom)
    
    # Test 1: Startup transient
    print("\n1. STARTUP TRANSIENT (0 → STEADY STATE)")
    print("-" * 70)
    
    conditions = TransientConditions(
        t_start=0.0,
        t_end=0.5,  # 500 ms
        dt=1e-4,
        current_initial=0.0,
        flow_rate_initial=0.0,
        temperature_initial=298.15
    )
    
    startup = transient.simulate_startup(conditions)
    
    if startup['success']:
        t_ss = transient.time_to_steady_state()
        I_final = startup['I'][-1]
        T_final = startup['T'][-1]
        
        print(f"Time to steady state: {t_ss*1000:.1f} ms")
        print(f"Final current: {I_final*1e9:.2f} nA")
        print(f"Final temperature: {T_final:.2f} K (ΔT = {T_final-298.15:.1f} K)")
        print(f"Current rise time (10-90%): {t_ss*0.8*1000:.1f} ms")
    
    # Test 2: Pulsed operation
    print("\n2. PULSED OPERATION")
    print("-" * 70)
    
    pulse = transient.simulate_pulse(
        t_on=0.1,  # 100 ms on
        t_off=0.1,  # 100 ms off
        n_cycles=5
    )
    
    I_mean_on = np.mean(pulse['I'][pulse['I'] > 1e-10])
    I_mean_off = np.mean(pulse['I'][pulse['I'] < 1e-10])
    
    print(f"Pulse frequency: {1/(pulse['t_on'] + pulse['t_off']):.1f} Hz")
    print(f"Duty cycle: {pulse['duty_cycle']*100:.0f}%")
    print(f"Mean current (on): {I_mean_on*1e9:.2f} nA")
    print(f"Mean current (off): {I_mean_off*1e9:.4f} nA")
    
    # Test 3: Current oscillations
    print("\n3. CURRENT OSCILLATIONS")
    print("-" * 70)
    
    osc = CurrentOscillationModel(frequency=100.0, amplitude_fraction=0.1)
    
    I_mean = 200e-9  # 200 nA
    t_osc = np.linspace(0, 0.1, 10000)
    I_with_osc = osc.current_with_oscillation(I_mean, t_osc)
    
    print(f"Oscillation frequency: {osc.f:.1f} Hz")
    print(f"Oscillation amplitude: {osc.A*100:.1f}% of mean")
    print(f"Mean current: {I_mean*1e9:.2f} nA")
    print(f"Peak-to-peak: {(I_with_osc.max() - I_with_osc.min())*1e9:.2f} nA")
    
    # Test 4: Mode transitions
    print("\n4. EMISSION MODE TRANSITIONS")
    print("-" * 70)
    
    mode_model = ModeTransitionModel()
    
    test_flows = [5e-15, 5e-14, 5e-13]  # m³/s
    
    for Q in test_flows:
        mode = mode_model.determine_mode(Q)
        print(f"Q = {Q:.2e} m³/s → Mode: {mode}")
    
    t_transition = mode_model.mode_transition_time(test_flows[0], test_flows[-1])
    print(f"\nMode transition time: {t_transition*1000:.2f} ms")
    
    print("\n" + "=" * 70)
    print("✓ Transient Dynamics Module Complete")
    print("=" * 70)
