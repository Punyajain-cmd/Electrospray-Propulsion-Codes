"""
Flow System Physics Module
===========================

Models the complete propellant feed system from tank to emitter.

Critical for hardware design - cannot design feed system without this.

Physics included:
1. Capillary flow resistance (Hagen-Poiseuille)
2. Hydraulic impedance
3. Pressure-flow relationships
4. Flow stability and oscillations
5. Tank pressurization
6. Flow rate control

References:
- Legge & Lozano (2011) - "Electrospray propulsion based on emitters microfabricated in porous metals"
- Demmons et al. (2018) - "Propellant feed system for ST7-DRS"
- Courtney et al. (2016) - "Ionic liquid ion source emitter arrays fabricated on bulk porous substrates"
"""

import numpy as np
from typing import Tuple, Optional, Dict
from dataclasses import dataclass


@dataclass
class FlowSystemGeometry:
    """Geometric parameters of flow system"""
    # Tank
    tank_volume: float  # m³
    tank_initial_pressure: float  # Pa (gauge)
    
    # Feed line
    capillary_length: float  # m
    capillary_diameter: float  # m
    
    # Emitter
    porous_length: float  # m (if using porous emitter)
    porous_permeability: Optional[float] = None  # m² (for porous media)
    
    # Optional: multiple capillaries in parallel
    num_parallel_lines: int = 1


class FlowSystemPhysics:
    """
    Complete flow system analysis
    
    Models pressure drop from tank through feed lines to emitter tip.
    Enables design of feed system for target flow rate.
    """
    
    def __init__(self, geometry: FlowSystemGeometry, liquid_props: Dict):
        self.geom = geometry
        self.liquid = liquid_props
        
        # Extract liquid properties
        self.mu = liquid_props['viscosity']  # Pa·s
        self.rho = liquid_props['density']  # kg/m³
        self.gamma = liquid_props['surface_tension']  # N/m
    
    def capillary_flow_resistance(self) -> float:
        """
        Hydraulic resistance of capillary feed line
        
        From Hagen-Poiseuille equation:
        R_h = 128μL/(πd⁴) for laminar flow in circular tube
        
        Returns:
            Hydraulic resistance (Pa·s/m³)
        """
        L = self.geom.capillary_length
        d = self.geom.capillary_diameter
        
        # Hagen-Poiseuille resistance
        R_h = 128 * self.mu * L / (np.pi * d**4)
        
        # Account for parallel lines
        if self.geom.num_parallel_lines > 1:
            R_h /= self.geom.num_parallel_lines
        
        return R_h
    
    def porous_media_resistance(self) -> float:
        """
        Hydraulic resistance through porous emitter substrate
        
        From Darcy's law:
        R_porous = μL/(κA)
        
        where κ is permeability, A is cross-sectional area
        
        Returns:
            Porous media resistance (Pa·s/m³)
        """
        if self.geom.porous_permeability is None:
            return 0.0  # Not a porous emitter
        
        L = self.geom.porous_length
        kappa = self.geom.porous_permeability
        A = np.pi * (self.geom.capillary_diameter / 2)**2
        
        return self.mu * L / (kappa * A)
    
    def total_hydraulic_resistance(self) -> float:
        """
        Total resistance: series combination
        
        R_total = R_capillary + R_porous
        """
        return self.capillary_flow_resistance() + self.porous_media_resistance()
    
    def pressure_drop(self, Q: float) -> float:
        """
        Pressure drop for given volumetric flow rate
        
        ΔP = R_h × Q
        
        Args:
            Q: Volumetric flow rate (m³/s)
        
        Returns:
            Pressure drop (Pa)
        """
        R_h = self.total_hydraulic_resistance()
        return R_h * Q
    
    def flow_rate_for_pressure(self, Delta_P: float) -> float:
        """
        Flow rate for given pressure drop
        
        Q = ΔP / R_h
        
        Args:
            Delta_P: Pressure drop (Pa)
        
        Returns:
            Volumetric flow rate (m³/s)
        """
        R_h = self.total_hydraulic_resistance()
        return Delta_P / R_h
    
    def reynolds_number(self, Q: float) -> float:
        """
        Reynolds number in feed capillary
        
        Re = ρvd/μ = 4ρQ/(πdμ)
        
        Re < 2300: Laminar (Hagen-Poiseuille valid)
        Re > 4000: Turbulent (need different model)
        """
        d = self.geom.capillary_diameter
        v = 4 * Q / (np.pi * d**2)  # Mean velocity
        
        return self.rho * v * d / self.mu
    
    def is_flow_laminar(self, Q: float) -> Tuple[bool, float]:
        """
        Check if flow is laminar
        
        Returns:
            (is_laminar, Re)
        """
        Re = self.reynolds_number(Q)
        return Re < 2300, Re
    
    def capillary_number(self, Q: float) -> float:
        """
        Capillary number: ratio of viscous to surface tension forces
        
        Ca = μv/γ
        
        Ca << 1: Surface tension dominates (stable meniscus)
        Ca >> 1: Viscous forces dominate (unstable)
        """
        d = self.geom.capillary_diameter
        v = 4 * Q / (np.pi * d**2)
        
        return self.mu * v / self.gamma
    
    def tank_blowdown_model(self, Q: float, t: float) -> Dict:
        """
        Model tank pressure decrease over time (blow-down mode)
        
        Assumes isothermal expansion of pressurant gas:
        P₁V₁ = P₂V₂
        
        Args:
            Q: Flow rate (m³/s) - assumed constant
            t: Time elapsed (s)
        
        Returns:
            Dict with P(t), V_liquid(t), V_gas(t)
        """
        V_tank = self.geom.tank_volume
        P_initial = self.geom.tank_initial_pressure
        
        # Volume of liquid consumed
        V_consumed = Q * t
        
        # Assuming initial tank half full
        V_liquid_initial = V_tank / 2
        V_gas_initial = V_tank / 2
        
        # Current volumes
        V_liquid = V_liquid_initial - V_consumed
        V_gas = V_tank - V_liquid
        
        # Check if tank empty
        if V_liquid <= 0:
            return {
                'pressure': 0,
                'V_liquid': 0,
                'V_gas': V_tank,
                'empty': True
            }
        
        # Isothermal gas expansion: P₁V₁ = P₂V₂
        P_current = P_initial * V_gas_initial / V_gas
        
        return {
            'pressure': P_current,
            'V_liquid': V_liquid,
            'V_gas': V_gas,
            'empty': False
        }
    
    def mission_duration(self, m_dot_required: float) -> Dict:
        """
        Calculate mission duration for given mass flow rate
        
        Args:
            m_dot_required: Required mass flow rate (kg/s)
        
        Returns:
            Dict with duration, propellant mass, etc.
        """
        # Total propellant mass (assume tank half full initially)
        V_propellant = self.geom.tank_volume / 2
        m_propellant = V_propellant * self.rho
        
        # Time to depletion
        t_mission = m_propellant / m_dot_required
        
        # Check if pressure adequate throughout
        Q_required = m_dot_required / self.rho
        Delta_P_required = self.pressure_drop(Q_required)
        
        # Final tank pressure (after propellant consumed)
        tank_final = self.tank_blowdown_model(Q_required, t_mission)
        
        return {
            'mission_duration': t_mission,
            'propellant_mass': m_propellant,
            'final_pressure': tank_final['pressure'],
            'pressure_adequate': tank_final['pressure'] > Delta_P_required,
            'margin': tank_final['pressure'] / Delta_P_required if Delta_P_required > 0 else np.inf
        }
    
    def design_feed_system(self, Q_target: float, P_tank_max: float = 200e3) -> Dict:
        """
        Design feed system for target flow rate
        
        Determines required capillary diameter and tank pressure.
        
        Args:
            Q_target: Target volumetric flow rate (m³/s)
            P_tank_max: Maximum allowable tank pressure (Pa)
        
        Returns:
            Design parameters
        """
        # Required pressure drop
        Delta_P_required = self.pressure_drop(Q_target)
        
        # Add safety margin (typically 2x)
        P_tank_required = Delta_P_required * 2.0
        
        # Check Reynolds number
        is_laminar, Re = self.is_flow_laminar(Q_target)
        
        # Check capillary number
        Ca = self.capillary_number(Q_target)
        
        # Assess design
        feasible = (P_tank_required < P_tank_max) and is_laminar and (Ca < 0.1)
        
        warnings = []
        if not is_laminar:
            warnings.append(f"Flow is turbulent (Re={Re:.0f}), model may be inaccurate")
        if Ca > 0.1:
            warnings.append(f"High capillary number (Ca={Ca:.3f}), meniscus may be unstable")
        if P_tank_required > P_tank_max:
            warnings.append(f"Required pressure ({P_tank_required/1e3:.1f} kPa) exceeds limit")
        
        return {
            'Q_target': Q_target,
            'Delta_P': Delta_P_required,
            'P_tank_required': P_tank_required,
            'P_tank_margin': P_tank_max / P_tank_required,
            'Re': Re,
            'Ca': Ca,
            'feasible': feasible,
            'warnings': warnings,
            'R_hydraulic': self.total_hydraulic_resistance()
        }


class FlowOscillationModel:
    """
    Models flow rate oscillations and instabilities
    
    Flow can oscillate due to:
    1. Meniscus dynamics
    2. Pressure fluctuations
    3. Electrohydrodynamic coupling
    """
    
    def __init__(self, flow_system: FlowSystemPhysics, emitter_capacitance: float):
        self.flow_sys = flow_system
        self.C_hydraulic = emitter_capacitance  # Hydraulic capacitance (m³/Pa)
    
    def natural_frequency(self) -> float:
        """
        Natural frequency of flow oscillations
        
        For hydraulic RC circuit:
        f = 1/(2π√(RC))
        
        Returns:
            Natural frequency (Hz)
        """
        R = self.flow_sys.total_hydraulic_resistance()
        C = self.C_hydraulic
        
        omega = 1 / np.sqrt(R * C)
        f = omega / (2 * np.pi)
        
        return f
    
    def damping_ratio(self, viscous_damping: float = 0.1) -> float:
        """
        Damping ratio of oscillations
        
        ζ < 1: Underdamped (oscillatory)
        ζ = 1: Critically damped
        ζ > 1: Overdamped
        
        Args:
            viscous_damping: Damping coefficient
        
        Returns:
            Damping ratio
        """
        # Simplified model
        return viscous_damping
    
    def transient_response(self, t_array: np.ndarray, 
                          Q_initial: float, Q_final: float) -> np.ndarray:
        """
        Flow rate response to step change
        
        Args:
            t_array: Time array (s)
            Q_initial: Initial flow rate (m³/s)
            Q_final: Final (steady-state) flow rate (m³/s)
        
        Returns:
            Q(t) flow rate array
        """
        f_n = self.natural_frequency()
        zeta = self.damping_ratio()
        omega_n = 2 * np.pi * f_n
        
        # Step response
        if zeta < 1:  # Underdamped
            omega_d = omega_n * np.sqrt(1 - zeta**2)
            Q_t = Q_final + (Q_initial - Q_final) * np.exp(-zeta * omega_n * t_array) * \
                  np.cos(omega_d * t_array)
        else:  # Overdamped or critically damped
            Q_t = Q_final + (Q_initial - Q_final) * np.exp(-omega_n * t_array)
        
        return Q_t


# ============================================================================
# INTEGRATION WITH EXISTING MODEL
# ============================================================================

class EnhancedEmitterWithFlow:
    """
    Wrapper that adds flow system to existing emitter model
    
    This preserves backward compatibility while adding new capability.
    """
    
    def __init__(self, base_emitter, flow_geometry: Optional[FlowSystemGeometry] = None):
        """
        Args:
            base_emitter: Existing SingleCapillaryEmitter instance
            flow_geometry: Optional flow system geometry
        """
        self.emitter = base_emitter
        
        if flow_geometry is not None:
            # Enable flow system modeling
            liquid_props = {
                'viscosity': base_emitter.liquid.viscosity,
                'density': base_emitter.liquid.density,
                'surface_tension': base_emitter.liquid.surface_tension,
            }
            self.flow_system = FlowSystemPhysics(flow_geometry, liquid_props)
            self.has_flow_system = True
        else:
            self.flow_system = None
            self.has_flow_system = False
    
    def calculate_with_flow_constraints(self):
        """
        Calculate emission with flow system constraints
        
        Returns:
            Dict with emission results AND flow system analysis
        """
        # Get base emission (existing model)
        self.emitter.calculate_emission()
        
        if not self.has_flow_system:
            # No flow system, return base results
            return {
                'current': self.emitter.I,
                'thrust': getattr(self.emitter, 'thrust', None),
                'flow_system_enabled': False
            }
        
        # With flow system: check if pressure is adequate
        Q_required = self.emitter.m_dot / self.emitter.liquid.density
        
        # Pressure drop in feed system
        Delta_P = self.flow_system.pressure_drop(Q_required)
        
        # Tank pressure
        tank_state = self.flow_system.tank_blowdown_model(Q_required, t=0)
        P_available = tank_state['pressure']
        
        # Check feasibility
        feasible = P_available > Delta_P
        
        return {
            'current': self.emitter.I,
            'thrust': getattr(self.emitter, 'thrust', None),
            'flow_system_enabled': True,
            'Q_required': Q_required,
            'P_available': P_available,
            'Delta_P': Delta_P,
            'feasible': feasible,
            'pressure_margin': P_available / Delta_P if Delta_P > 0 else np.inf,
            'Re': self.flow_system.reynolds_number(Q_required),
        }


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Flow System Physics Module - Demonstration")
    print("=" * 70)
    
    # Define flow system geometry
    flow_geom = FlowSystemGeometry(
        tank_volume=100e-6,  # 100 mL
        tank_initial_pressure=200e3,  # 200 kPa
        capillary_length=0.05,  # 5 cm feed line
        capillary_diameter=100e-6,  # 100 μm inner diameter
        porous_length=0.001,  # 1 mm porous substrate
        porous_permeability=1e-15,  # 1 μm² equivalent
    )
    
    # Liquid properties (EMI-Im example)
    liquid_props = {
        'viscosity': 0.034,  # Pa·s
        'density': 1520,  # kg/m³
        'surface_tension': 0.042,  # N/m
    }
    
    # Create flow system
    flow_sys = FlowSystemPhysics(flow_geom, liquid_props)
    
    # Test: Design for target flow rate
    print("\n1. FLOW SYSTEM DESIGN")
    print("-" * 70)
    
    m_dot_target = 1e-10  # 100 pg/s (typical for single emitter)
    Q_target = m_dot_target / liquid_props['density']
    
    design = flow_sys.design_feed_system(Q_target)
    
    print(f"Target flow rate: {m_dot_target*1e12:.2f} pg/s")
    print(f"Pressure drop: {design['Delta_P']/1e3:.2f} kPa")
    print(f"Required tank pressure: {design['P_tank_required']/1e3:.2f} kPa")
    print(f"Pressure margin: {design['P_tank_margin']:.2f}x")
    print(f"Reynolds number: {design['Re']:.2f} (laminar: <2300)")
    print(f"Capillary number: {design['Ca']:.4f} (stable: <0.1)")
    print(f"Feasible: {design['feasible']}")
    
    if design['warnings']:
        print("\nWarnings:")
        for w in design['warnings']:
            print(f"  ⚠ {w}")
    
    # Test: Mission duration
    print("\n2. MISSION DURATION ANALYSIS")
    print("-" * 70)
    
    mission = flow_sys.mission_duration(m_dot_target)
    
    print(f"Propellant mass: {mission['propellant_mass']*1e6:.2f} mg")
    print(f"Mission duration: {mission['mission_duration']/3600:.2f} hours")
    print(f"Final tank pressure: {mission['final_pressure']/1e3:.2f} kPa")
    print(f"Pressure adequate: {mission['pressure_adequate']}")
    print(f"Pressure margin: {mission['margin']:.2f}x")
    
    # Test: Hydraulic resistance
    print("\n3. HYDRAULIC RESISTANCE BREAKDOWN")
    print("-" * 70)
    
    R_cap = flow_sys.capillary_flow_resistance()
    R_porous = flow_sys.porous_media_resistance()
    R_total = flow_sys.total_hydraulic_resistance()
    
    print(f"Capillary resistance: {R_cap:.3e} Pa·s/m³")
    print(f"Porous media resistance: {R_porous:.3e} Pa·s/m³")
    print(f"Total resistance: {R_total:.3e} Pa·s/m³")
    print(f"Porous fraction: {R_porous/R_total*100:.1f}%")
    
    print("\n" + "=" * 70)
    print("✓ Flow System Module Complete")
    print("=" * 70)
