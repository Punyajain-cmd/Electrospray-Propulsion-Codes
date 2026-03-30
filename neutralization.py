"""
Beam Neutralization Module
===========================

Models neutralization of ion beam to prevent spacecraft charging.

CRITICAL for spacecraft operation - without neutralization, spacecraft
would charge to kilovolts and arc/damage electronics.

Physics included:
1. Electron emission (thermionic, field emission, hollow cathode)
2. Beam-plasma coupling
3. Spacecraft charging analysis
4. Neutralizer design
5. Plume neutrality

References:
- Marrese et al. (2002) - "Electric Propulsion Plume Characterization"
- Berg & Rovey (2015) - "Electrospray thruster performance and lifetime"
- Reza et al. (2020) - "Bipolar colloid thrusters"
"""

import numpy as np
from typing import Dict, Tuple, Optional
from dataclasses import dataclass


@dataclass
class NeutralizerGeometry:
    """Geometric parameters of neutralizer"""
    # Cathode
    cathode_diameter: float  # m
    cathode_length: float  # m
    cathode_material: str  # 'tungsten', 'lanthanum_hexaboride', etc.
    
    # Emission surface
    emission_area: float  # m²
    work_function: float  # eV
    
    # Position relative to thruster
    distance_from_thruster: float  # m


class ElectronEmissionModel:
    """
    Models electron emission from neutralizer cathode
    
    Mechanisms:
    1. Thermionic emission (Richardson-Dushman)
    2. Field-enhanced thermionic (Schottky)
    3. Field emission (Fowler-Nordheim)
    """
    
    def __init__(self, work_function: float, emission_area: float,
                 material: str = "tungsten"):
        """
        Args:
            work_function: Material work function (eV)
            emission_area: Effective emission area (m²)
            material: Cathode material
        """
        self.phi = work_function  # eV
        self.A_em = emission_area  # m²
        self.material = material
        
        # Richardson constant
        self.A_0 = 1.2e6  # A/(m²·K²) - universal
        
        # Physical constants
        self.e = 1.602e-19  # C
        self.k_B = 1.381e-23  # J/K
        self.m_e = 9.109e-31  # kg
        self.h = 6.626e-34  # J·s
    
    def thermionic_current(self, T: float) -> float:
        """
        Thermionic emission current (Richardson-Dushman equation)
        
        J = A₀T² exp(-eφ/kT)
        I = J × A
        
        Args:
            T: Cathode temperature (K)
        
        Returns:
            Emission current (A)
        """
        # Current density
        J = self.A_0 * T**2 * np.exp(-self.e * self.phi / (self.k_B * T))
        
        # Total current
        I = J * self.A_em
        
        return I
    
    def schottky_enhancement(self, E_field: float, T: float) -> float:
        """
        Field-enhanced thermionic emission (Schottky effect)
        
        Effective work function reduced by electric field:
        φ_eff = φ - √(e³E/(4πε₀))
        
        Args:
            E_field: Electric field at surface (V/m)
            T: Temperature (K)
        
        Returns:
            Enhanced emission current (A)
        """
        epsilon_0 = 8.854e-12  # F/m
        
        # Work function reduction (eV)
        Delta_phi = np.sqrt(self.e**3 * E_field / (4 * np.pi * epsilon_0)) / self.e
        
        # Effective work function
        phi_eff = self.phi - Delta_phi
        
        # Enhanced emission
        J = self.A_0 * T**2 * np.exp(-self.e * phi_eff / (self.k_B * T))
        I = J * self.A_em
        
        return I
    
    def field_emission_current(self, E_field: float) -> float:
        """
        Pure field emission (Fowler-Nordheim)
        
        J = (e³E²)/(16π²ħφ) exp(-4√(2m)φ^(3/2)/(3ħeE))
        
        Dominant at high fields (>10⁹ V/m), low temperatures
        
        Args:
            E_field: Electric field (V/m)
        
        Returns:
            Field emission current (A)
        """
        # Simplified Fowler-Nordheim
        # J ≈ a E² exp(-b/E)  where a, b are material-dependent
        
        a = 1.54e-6 / self.phi  # A·m²/V²
        b = 6.83e9 * self.phi**(3/2)  # V/m
        
        J = a * E_field**2 * np.exp(-b / E_field)
        I = J * self.A_em
        
        return I
    
    def total_emission_current(self, T: float, E_field: float) -> Dict:
        """
        Total emission from all mechanisms
        
        Returns:
            Dict with breakdown by mechanism
        """
        I_thermionic = self.thermionic_current(T)
        I_schottky = self.schottky_enhancement(E_field, T) - I_thermionic
        I_field = self.field_emission_current(E_field)
        
        # At high fields/temps, mechanisms compete
        # Use max of (thermionic, schottky, field)
        I_total = max(I_thermionic, I_thermionic + I_schottky, I_field)
        
        return {
            'total': I_total,
            'thermionic': I_thermionic,
            'schottky_enhancement': I_schottky,
            'field_emission': I_field,
            'dominant_mechanism': 'thermionic' if I_thermionic > I_field else 'field'
        }


class SpacecraftChargingModel:
    """
    Models spacecraft charging in absence of neutralization
    
    Critical for assessing neutralizer requirement.
    """
    
    def __init__(self, spacecraft_area: float, spacecraft_capacitance: float):
        """
        Args:
            spacecraft_area: Total surface area (m²)
            spacecraft_capacitance: Self-capacitance (F)
        """
        self.A_sc = spacecraft_area
        self.C_sc = spacecraft_capacitance
        
        # Physical constants
        self.e = 1.602e-19  # C
        self.epsilon_0 = 8.854e-12  # F/m
    
    def charging_rate(self, I_beam: float, I_return: float = 0.0) -> float:
        """
        Rate of spacecraft potential increase
        
        dV/dt = (I_beam - I_return) / C_sc
        
        Args:
            I_beam: Ion beam current (A) - leaving spacecraft
            I_return: Return current (A) - electrons collected
        
        Returns:
            dV/dt (V/s)
        """
        I_net = I_beam - I_return
        return I_net / self.C_sc
    
    def time_to_voltage(self, I_beam: float, V_target: float, 
                        I_return: float = 0.0) -> float:
        """
        Time to charge to target voltage
        
        Args:
            I_beam: Beam current (A)
            V_target: Target voltage (V)
            I_return: Return current (A)
        
        Returns:
            Time (s)
        """
        dV_dt = self.charging_rate(I_beam, I_return)
        
        if dV_dt <= 0:
            return np.inf  # Never reaches (neutral or discharging)
        
        return V_target / dV_dt
    
    def equilibrium_potential(self, I_beam: float, 
                              ambient_ne: float, ambient_Te: float) -> float:
        """
        Equilibrium floating potential in plasma environment
        
        Balance: I_beam = I_electron_collection
        
        Args:
            I_beam: Beam current (A)
            ambient_ne: Ambient electron density (m⁻³)
            ambient_Te: Ambient electron temperature (eV)
        
        Returns:
            Equilibrium potential (V)
        """
        # Electron thermal velocity
        m_e = 9.109e-31  # kg
        v_th = np.sqrt(2 * self.e * ambient_Te / m_e)
        
        # Random electron current to surface
        I_e_thermal = 0.25 * self.e * ambient_ne * v_th * self.A_sc
        
        # If no neutralizer, spacecraft charges positive until
        # electron collection = ion emission
        
        # Sheath potential to collect I_beam worth of electrons
        # I_e = I_e_thermal * exp(eV/kTe)
        # V = (kTe/e) ln(I_beam / I_e_thermal)
        
        if I_beam <= I_e_thermal:
            return 0.0  # No significant charging
        
        V_float = ambient_Te * np.log(I_beam / I_e_thermal)
        
        return V_float


class NeutralizerDesignModel:
    """
    Design neutralizer to match thruster current
    """
    
    def __init__(self, thruster_current: float):
        self.I_thruster = thruster_current
    
    def required_electron_current(self, margin: float = 1.2) -> float:
        """
        Required electron emission current
        
        Should exceed thruster current by margin for stability
        
        Args:
            margin: Safety margin (e.g., 1.2 = 20% over-emission)
        
        Returns:
            Required I_electron (A)
        """
        return self.I_thruster * margin
    
    def design_thermionic_cathode(self, T_max: float = 1800) -> Dict:
        """
        Design thermionic cathode for required current
        
        Args:
            T_max: Maximum allowable cathode temperature (K)
        
        Returns:
            Design parameters
        """
        I_required = self.required_electron_current()
        
        # Work functions for common materials
        materials = {
            'tungsten': 4.5,  # eV
            'tantalum': 4.25,
            'lanthanum_hexaboride': 2.7,
            'barium_oxide': 1.5,
        }
        
        designs = {}
        
        for material, phi in materials.items():
            # Required emission current density at T_max
            A_0 = 1.2e6
            e = 1.602e-19
            k_B = 1.381e-23
            
            J_max = A_0 * T_max**2 * np.exp(-e * phi / (k_B * T_max))
            
            # Required area
            A_required = I_required / J_max
            
            # Equivalent radius (assuming cylindrical)
            r_cathode = np.sqrt(A_required / (2 * np.pi * 0.01))  # L = 1 cm assumed
            
            designs[material] = {
                'work_function': phi,
                'current_density': J_max,
                'required_area': A_required,
                'cathode_radius': r_cathode,
                'temperature': T_max,
                'feasible': A_required < 1e-4  # < 1 cm²
            }
        
        return designs
    
    def power_requirement(self, T_cathode: float, A_cathode: float,
                         emissivity: float = 0.3) -> float:
        """
        Heater power requirement for thermionic cathode
        
        P = σεAT⁴ (radiative cooling)
        
        Args:
            T_cathode: Operating temperature (K)
            A_cathode: Cathode surface area (m²)
            emissivity: Surface emissivity
        
        Returns:
            Power (W)
        """
        sigma = 5.67e-8  # W/(m²·K⁴) - Stefan-Boltzmann constant
        
        P_rad = sigma * emissivity * A_cathode * T_cathode**4
        
        # Add efficiency factor (heater not 100% efficient)
        eta_heater = 0.5
        
        P_required = P_rad / eta_heater
        
        return P_required


class BipolarEmissionModel:
    """
    Model for bipolar emission (alternating positive/negative ions)
    
    Eliminates need for separate neutralizer - thruster emits both
    polarities in sequence, achieving net neutrality.
    """
    
    def __init__(self, I_positive: float, I_negative: float, 
                 f_alternation: float):
        """
        Args:
            I_positive: Positive ion current (A)
            I_negative: Negative ion current (A)
            f_alternation: Alternation frequency (Hz)
        """
        self.I_pos = I_positive
        self.I_neg = I_negative
        self.f_alt = f_alternation
    
    def charge_balance_per_cycle(self) -> float:
        """
        Net charge emitted per alternation cycle
        
        Q_net = (I_pos - I_neg) / f
        
        Returns:
            Net charge per cycle (C)
        """
        Q_cycle = (self.I_pos - self.I_neg) / self.f_alt
        return Q_cycle
    
    def is_neutral(self, tolerance: float = 0.05) -> bool:
        """
        Check if emission is approximately neutral
        
        Args:
            tolerance: Allowable imbalance fraction (0.05 = 5%)
        
        Returns:
            True if balanced within tolerance
        """
        I_avg = (self.I_pos + self.I_neg) / 2
        imbalance = abs(self.I_pos - self.I_neg) / I_avg
        
        return imbalance < tolerance


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Beam Neutralization Module - Demonstration")
    print("=" * 70)
    
    # Test 1: Electron emission
    print("\n1. ELECTRON EMISSION")
    print("-" * 70)
    
    emitter = ElectronEmissionModel(
        work_function=2.7,  # eV - LaB₆
        emission_area=1e-6,  # 1 mm²
        material="lanthanum_hexaboride"
    )
    
    T_cathode = 1500  # K
    E_field = 1e7  # 10 MV/m
    
    emission = emitter.total_emission_current(T_cathode, E_field)
    
    print(f"Cathode temperature: {T_cathode} K")
    print(f"Work function: {emitter.phi} eV")
    print(f"Emission area: {emitter.A_em*1e6:.2f} mm²")
    print(f"\nEmission currents:")
    print(f"  Thermionic: {emission['thermionic']*1e6:.2f} μA")
    print(f"  Schottky enhancement: {emission['schottky_enhancement']*1e6:.2f} μA")
    print(f"  Field emission: {emission['field_emission']*1e6:.2f} μA")
    print(f"  TOTAL: {emission['total']*1e6:.2f} μA")
    print(f"  Dominant: {emission['dominant_mechanism']}")
    
    # Test 2: Spacecraft charging
    print("\n2. SPACECRAFT CHARGING (Without Neutralizer)")
    print("-" * 70)
    
    charging = SpacecraftChargingModel(
        spacecraft_area=1.0,  # 1 m² surface
        spacecraft_capacitance=100e-12  # 100 pF
    )
    
    I_thruster = 200e-9  # 200 nA
    
    dV_dt = charging.charging_rate(I_thruster, I_return=0)
    t_1kV = charging.time_to_voltage(I_thruster, 1000, I_return=0)
    
    print(f"Thruster current: {I_thruster*1e9:.0f} nA")
    print(f"Spacecraft capacitance: {charging.C_sc*1e12:.0f} pF")
    print(f"Charging rate: {dV_dt:.2f} V/s")
    print(f"Time to 1 kV: {t_1kV:.1f} seconds")
    print(f"⚠ Neutralizer REQUIRED to prevent charging!")
    
    # Test 3: Neutralizer design
    print("\n3. NEUTRALIZER DESIGN")
    print("-" * 70)
    
    designer = NeutralizerDesignModel(thruster_current=200e-9)
    
    I_e_required = designer.required_electron_current()
    print(f"Required electron current: {I_e_required*1e9:.0f} nA (with 20% margin)")
    
    designs = designer.design_thermionic_cathode(T_max=1800)
    
    print(f"\nCathode material options:")
    for material, design in designs.items():
        print(f"\n  {material}:")
        print(f"    Work function: {design['work_function']:.2f} eV")
        print(f"    Current density: {design['current_density']:.2e} A/m²")
        print(f"    Required area: {design['required_area']*1e6:.2f} mm²")
        print(f"    Cathode radius: {design['cathode_radius']*1e3:.2f} mm")
        print(f"    Feasible: {design['feasible']}")
    
    # Select LaB₆ (good balance)
    print(f"\n  RECOMMENDED: Lanthanum hexaboride (LaB₆)")
    A_lab6 = designs['lanthanum_hexaboride']['required_area']
    P_heater = designer.power_requirement(1800, A_lab6)
    print(f"  Heater power: {P_heater*1e3:.2f} mW")
    
    # Test 4: Bipolar operation (alternative)
    print("\n4. BIPOLAR OPERATION (Alternative to Neutralizer)")
    print("-" * 70)
    
    bipolar = BipolarEmissionModel(
        I_positive=200e-9,
        I_negative=195e-9,  # Slight imbalance
        f_alternation=100  # 100 Hz
    )
    
    Q_net = bipolar.charge_balance_per_cycle()
    is_neutral = bipolar.is_neutral(tolerance=0.05)
    
    print(f"Positive current: {bipolar.I_pos*1e9:.0f} nA")
    print(f"Negative current: {bipolar.I_neg*1e9:.0f} nA")
    print(f"Alternation frequency: {bipolar.f_alt:.0f} Hz")
    print(f"Net charge/cycle: {Q_net*1e12:.2f} pC")
    print(f"Neutral (within 5%): {is_neutral}")
    
    print("\n" + "=" * 70)
    print("✓ Neutralization Module Complete")
    print("=" * 70)
