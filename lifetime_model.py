"""
Lifetime and Degradation Module
================================

Models long-term performance degradation and predicts mission lifetime.

Critical for mission planning - need to know how long thruster will last.

Physics included:
1. Propellant depletion tracking
2. Emitter erosion/degradation
3. Wetting property changes
4. Contamination buildup
5. Performance drift over time

References:
- Freeman et al. (2011) - "Lifetime testing of porous emitter microthruster"
- Petro et al. (2022) - "Lifetime testing and degradation mechanisms"
- Courtney et al. (2016) - "Long duration testing of ionic liquid ion sources"
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import warnings


@dataclass
class MissionProfile:
    """Mission operating profile"""
    total_mission_time: float  # s - total mission duration
    duty_cycle: float  # 0-1 - fraction of time thruster is on
    average_current: float  # A - time-averaged current
    average_thrust: float  # N - time-averaged thrust
    num_on_off_cycles: int  # Number of startup/shutdown cycles


class PropellantDepletionModel:
    """
    Track propellant consumption and predict depletion time
    """
    
    def __init__(self, tank_capacity: float, liquid_density: float):
        """
        Args:
            tank_capacity: Total propellant volume (m³)
            liquid_density: Propellant density (kg/m³)
        """
        self.V_total = tank_capacity
        self.rho = liquid_density
        self.m_total = tank_capacity * liquid_density
        
        # Track consumption history
        self.t_history = [0.0]
        self.m_consumed_history = [0.0]
        self.m_remaining_history = [self.m_total]
    
    def mass_remaining(self, m_dot: float, t_elapsed: float) -> float:
        """
        Calculate remaining propellant mass
        
        Args:
            m_dot: Mass flow rate (kg/s)
            t_elapsed: Time elapsed (s)
        
        Returns:
            Mass remaining (kg)
        """
        m_consumed = m_dot * t_elapsed
        m_remaining = self.m_total - m_consumed
        
        return max(m_remaining, 0.0)
    
    def depletion_time(self, m_dot: float) -> float:
        """
        Time until propellant depleted
        
        Args:
            m_dot: Constant mass flow rate (kg/s)
        
        Returns:
            Time to depletion (s)
        """
        if m_dot <= 0:
            return np.inf
        
        return self.m_total / m_dot
    
    def depletion_time_with_duty_cycle(self, m_dot_on: float, duty_cycle: float) -> float:
        """
        Time until depletion with pulsed operation
        
        Args:
            m_dot_on: Flow rate when thruster is on (kg/s)
            duty_cycle: Fraction of time thruster is on (0-1)
        
        Returns:
            Mission duration (s)
        """
        m_dot_avg = m_dot_on * duty_cycle
        return self.depletion_time(m_dot_avg)
    
    def add_consumption_record(self, t: float, m_consumed: float):
        """Record consumption at time t"""
        self.t_history.append(t)
        self.m_consumed_history.append(m_consumed)
        self.m_remaining_history.append(self.m_total - m_consumed)
    
    def get_consumption_rate(self) -> float:
        """
        Estimate current consumption rate from history
        
        Returns:
            dm/dt (kg/s)
        """
        if len(self.t_history) < 2:
            return 0.0
        
        # Linear fit to recent data
        recent_n = min(10, len(self.t_history))
        t_recent = np.array(self.t_history[-recent_n:])
        m_recent = np.array(self.m_consumed_history[-recent_n:])
        
        # Fit: m = m_dot * t
        m_dot = np.polyfit(t_recent, m_recent, 1)[0]
        
        return m_dot


class EmitterDegradationModel:
    """
    Models physical degradation of emitter structure
    
    Mechanisms:
    1. Sputtering erosion (ion bombardment)
    2. Chemical etching
    3. Mechanical stress (thermal cycling)
    4. Surface reconstruction
    """
    
    def __init__(self, emitter_material: str = "porous_tungsten"):
        self.material = emitter_material
        
        # Degradation rate constants (empirical, from lifetime tests)
        self.degradation_rates = {
            'porous_tungsten': 1e-18,  # m³/C - volume eroded per coulomb
            'stainless_steel': 5e-18,
            'silicon': 3e-18,
            'nickel': 4e-18,
        }
        
        self.k_erosion = self.degradation_rates.get(emitter_material, 2e-18)
    
    def erosion_depth(self, I: float, t: float) -> float:
        """
        Erosion depth from ion bombardment
        
        Depth ≈ k × Q  where Q = ∫I dt (total charge)
        
        Args:
            I: Current (A)
            t: Time (s)
        
        Returns:
            Erosion depth (m)
        """
        Q_total = I * t  # Total charge (C)
        
        # Volume eroded
        V_eroded = self.k_erosion * Q_total
        
        # Assume uniform erosion over area ~πr²
        r_tip = 15e-6  # m - typical tip radius
        A_tip = np.pi * r_tip**2
        
        depth = V_eroded / A_tip
        
        return depth
    
    def tip_radius_change(self, I: float, t: float, r_initial: float) -> float:
        """
        Change in tip radius due to erosion
        
        Args:
            I: Current (A)
            t: Time (s)
            r_initial: Initial tip radius (m)
        
        Returns:
            New tip radius (m)
        """
        depth = self.erosion_depth(I, t)
        
        # Erosion increases tip radius (blunting)
        r_new = r_initial + depth
        
        return r_new
    
    def performance_degradation_factor(self, I: float, t: float) -> float:
        """
        Factor by which performance degrades (0-1)
        
        1.0 = no degradation
        0.0 = complete failure
        
        Empirical model based on tip radius change
        """
        r_initial = 15e-6
        r_current = self.tip_radius_change(I, t, r_initial)
        
        # Performance degrades as tip blunts
        # Threshold: 2x radius increase = 50% performance
        degradation_factor = 1.0 / (1.0 + (r_current / r_initial - 1.0))
        
        return np.clip(degradation_factor, 0.0, 1.0)
    
    def estimated_lifetime(self, I: float, threshold: float = 0.7) -> float:
        """
        Estimate emitter lifetime
        
        Args:
            I: Operating current (A)
            threshold: Performance threshold (0-1) below which emitter fails
        
        Returns:
            Lifetime (s)
        """
        # Binary search for time when degradation factor = threshold
        t_min, t_max = 0, 1e10  # s
        
        for _ in range(50):  # iterations
            t_mid = (t_min + t_max) / 2
            factor = self.performance_degradation_factor(I, t_mid)
            
            if factor > threshold:
                t_min = t_mid
            else:
                t_max = t_mid
        
        return t_mid


class WettingDegradationModel:
    """
    Models degradation of wetting properties
    
    Over time, emitter surface can become contaminated or chemically
    altered, affecting wetting and meniscus stability.
    """
    
    def __init__(self, initial_contact_angle: float = 20.0):
        """
        Args:
            initial_contact_angle: Initial contact angle (degrees)
        """
        self.theta_0 = initial_contact_angle  # degrees
        
        # Degradation time constant (empirical)
        self.tau_wetting = 1e7  # s (~115 days)
    
    def contact_angle(self, t: float) -> float:
        """
        Contact angle evolution with time
        
        θ(t) = θ₀ + Δθ(1 - exp(-t/τ))
        
        Wetting degrades (angle increases) over time
        
        Args:
            t: Time (s)
        
        Returns:
            Contact angle (degrees)
        """
        Delta_theta = 30.0  # Saturation increase (degrees)
        
        theta_t = self.theta_0 + Delta_theta * (1 - np.exp(-t / self.tau_wetting))
        
        return theta_t
    
    def is_wetting_adequate(self, t: float) -> bool:
        """
        Check if wetting is still adequate for operation
        
        Threshold: θ < 60° for stable wetting
        """
        theta = self.contact_angle(t)
        return theta < 60.0
    
    def wetting_lifetime(self) -> float:
        """
        Time until wetting failure (θ > 60°)
        
        Returns:
            Lifetime (s)
        """
        # Solve: θ₀ + Δθ(1 - exp(-t/τ)) = 60
        # exp(-t/τ) = 1 - (60 - θ₀)/Δθ
        # t = -τ ln(1 - (60 - θ₀)/Δθ)
        
        Delta_theta = 30.0
        
        if self.theta_0 >= 60:
            return 0.0  # Already failed
        
        fraction = (60 - self.theta_0) / Delta_theta
        
        if fraction >= 1.0:
            return np.inf  # Never fails
        
        t_fail = -self.tau_wetting * np.log(1 - fraction)
        
        return t_fail


class ComprehensiveLifetimeModel:
    """
    Combines all degradation mechanisms to predict overall lifetime
    """
    
    def __init__(self, propellant_capacity: float, liquid_density: float,
                 emitter_material: str = "porous_tungsten"):
        self.propellant_model = PropellantDepletionModel(propellant_capacity, liquid_density)
        self.erosion_model = EmitterDegradationModel(emitter_material)
        self.wetting_model = WettingDegradationModel()
    
    def predict_lifetime(self, mission: MissionProfile) -> Dict:
        """
        Predict mission lifetime considering all degradation modes
        
        Args:
            mission: Mission operating profile
        
        Returns:
            Dict with lifetime predictions and limiting factor
        """
        # Average mass flow rate
        I_avg = mission.average_current
        psi = 2.48
        gamma = 0.042
        K = 1.5
        rho = 1520
        
        # Estimate average mass flow from current
        # I = psi√(γKm/ρ) → m = (I/psi)²(ρ/γK)
        m_dot_avg = (I_avg / psi)**2 * (rho / (gamma * K))
        
        # 1. Propellant depletion lifetime
        t_propellant = self.propellant_model.depletion_time_with_duty_cycle(
            m_dot_avg / mission.duty_cycle,  # Flow when on
            mission.duty_cycle
        )
        
        # 2. Erosion lifetime
        t_erosion = self.erosion_model.estimated_lifetime(I_avg, threshold=0.7)
        
        # 3. Wetting lifetime
        t_wetting = self.wetting_model.wetting_lifetime()
        
        # Overall lifetime = minimum of all
        lifetimes = {
            'propellant': t_propellant,
            'erosion': t_erosion,
            'wetting': t_wetting,
        }
        
        t_total = min(lifetimes.values())
        limiting_factor = min(lifetimes, key=lifetimes.get)
        
        return {
            'total_lifetime': t_total,
            'limiting_factor': limiting_factor,
            'propellant_lifetime': t_propellant,
            'erosion_lifetime': t_erosion,
            'wetting_lifetime': t_wetting,
            'lifetimes_hours': {k: v/3600 for k, v in lifetimes.items()},
            'mission_achievable': t_total >= mission.total_mission_time,
            'margin': t_total / mission.total_mission_time if mission.total_mission_time > 0 else np.inf
        }


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Lifetime and Degradation Module - Demonstration")
    print("=" * 70)
    
    # Define mission
    mission = MissionProfile(
        total_mission_time=365 * 24 * 3600,  # 1 year
        duty_cycle=0.5,  # On 50% of the time
        average_current=200e-9,  # 200 nA average
        average_thrust=1e-6,  # 1 μN average
        num_on_off_cycles=10000
    )
    
    # Create lifetime model
    lifetime_model = ComprehensiveLifetimeModel(
        propellant_capacity=100e-6,  # 100 mL
        liquid_density=1520,  # kg/m³
        emitter_material="porous_tungsten"
    )
    
    # Predict lifetime
    print("\n1. LIFETIME PREDICTION")
    print("-" * 70)
    
    prediction = lifetime_model.predict_lifetime(mission)
    
    print(f"Mission duration: {mission.total_mission_time/3600/24:.0f} days")
    print(f"Duty cycle: {mission.duty_cycle*100:.0f}%")
    print(f"Average current: {mission.average_current*1e9:.0f} nA")
    
    print(f"\nLifetime Predictions:")
    print(f"  Propellant: {prediction['propellant_lifetime']/3600/24:.1f} days")
    print(f"  Erosion: {prediction['erosion_lifetime']/3600/24:.1f} days")
    print(f"  Wetting: {prediction['wetting_lifetime']/3600/24:.1f} days")
    
    print(f"\nTotal lifetime: {prediction['total_lifetime']/3600/24:.1f} days")
    print(f"Limiting factor: {prediction['limiting_factor']}")
    print(f"Mission achievable: {prediction['mission_achievable']}")
    print(f"Lifetime margin: {prediction['margin']:.2f}x")
    
    # Test 2: Propellant tracking
    print("\n2. PROPELLANT CONSUMPTION")
    print("-" * 70)
    
    m_dot = 1e-10  # kg/s
    
    for t_days in [0, 30, 90, 180, 365]:
        t_sec = t_days * 24 * 3600
        m_remaining = lifetime_model.propellant_model.mass_remaining(m_dot, t_sec)
        percent_remaining = m_remaining / lifetime_model.propellant_model.m_total * 100
        
        print(f"Day {t_days:3d}: {m_remaining*1e6:.2f} mg remaining ({percent_remaining:.1f}%)")
    
    # Test 3: Erosion
    print("\n3. EMITTER EROSION")
    print("-" * 70)
    
    I = 200e-9  # A
    
    for t_days in [0, 100, 500, 1000]:
        t_sec = t_days * 24 * 3600
        depth = lifetime_model.erosion_model.erosion_depth(I, t_sec)
        r_init = 15e-6
        r_new = lifetime_model.erosion_model.tip_radius_change(I, t_sec, r_init)
        factor = lifetime_model.erosion_model.performance_degradation_factor(I, t_sec)
        
        print(f"Day {t_days:4d}: Erosion = {depth*1e9:.2f} nm, "
              f"Tip radius = {r_new*1e6:.2f} μm, "
              f"Performance = {factor*100:.1f}%")
    
    # Test 4: Wetting degradation
    print("\n4. WETTING DEGRADATION")
    print("-" * 70)
    
    for t_days in [0, 50, 100, 200]:
        t_sec = t_days * 24 * 3600
        theta = lifetime_model.wetting_model.contact_angle(t_sec)
        adequate = lifetime_model.wetting_model.is_wetting_adequate(t_sec)
        
        print(f"Day {t_days:3d}: Contact angle = {theta:.1f}°, "
              f"Adequate: {adequate}")
    
    t_wetting_fail = lifetime_model.wetting_model.wetting_lifetime()
    print(f"\nWetting failure time: {t_wetting_fail/3600/24:.0f} days")
    
    print("\n" + "=" * 70)
    print("✓ Lifetime Module Complete")
    print("=" * 70)
