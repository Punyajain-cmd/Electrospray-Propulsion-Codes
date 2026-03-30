"""
Comprehensive Electrospray Thruster Simulation Model
Based on physics from 4 key research papers

Features:
- Single and multi-capillary configurations
- Complete physics including self-heating
- Ion cluster formation modeling
- Performance optimization
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import Tuple, Dict, List
import warnings
warnings.filterwarnings('ignore')

# Physical Constants
EPSILON_0 = 8.854187817e-12  # F/m
E_CHARGE = 1.602176634e-19   # C
K_B = 1.380649e-23           # J/K
G_0 = 9.81                   # m/s²
AMU = 1.66053906660e-27      # kg

@dataclass
class IonicLiquid:
    """Properties of ionic liquid propellants"""
    name: str
    conductivity: float        # S/m
    surface_tension: float     # N/m
    viscosity: float           # Pa·s
    density: float             # kg/m³
    permittivity: float        # relative
    M_cation: float            # kg/mol
    M_anion: float             # kg/mol
    
    # Temperature dependence (Arrhenius-type)
    K_0: float = None          # Pre-exponential factor for conductivity
    E_a_K: float = None        # Activation energy for conductivity (J/mol)
    mu_0: float = None         # Pre-exponential factor for viscosity
    E_a_mu: float = None       # Activation energy for viscosity (J/mol)
    
    def __post_init__(self):
        """Set default temperature coefficients if not provided"""
        R_gas = 8.314  # J/(mol·K)
        T_ref = 298.15  # K
        
        if self.K_0 is None:
            # Typical activation energy ~0.3 eV for ionic liquids
            self.E_a_K = 0.3 * E_CHARGE * 6.022e23  # J/mol
            self.K_0 = self.conductivity * np.exp(self.E_a_K / (R_gas * T_ref))
        
        if self.mu_0 is None:
            # Typical activation energy ~0.4 eV for viscosity
            self.E_a_mu = 0.4 * E_CHARGE * 6.022e23  # J/mol
            self.mu_0 = self.viscosity * np.exp(-self.E_a_mu / (R_gas * T_ref))
    
    def conductivity_T(self, T: float) -> float:
        """Temperature-dependent conductivity"""
        R_gas = 8.314
        return self.K_0 * np.exp(-self.E_a_K / (R_gas * T))
    
    def viscosity_T(self, T: float) -> float:
        """Temperature-dependent viscosity"""
        R_gas = 8.314
        return self.mu_0 * np.exp(self.E_a_mu / (R_gas * T))

# Define common ionic liquids
EMI_IM = IonicLiquid(
    name="EMI-Im",
    conductivity=0.92,
    surface_tension=0.036,
    viscosity=0.0396,
    density=1520.0,
    permittivity=13.8,
    M_cation=111.2e-3,
    M_anion=280.1e-3
)

EMI_BF4 = IonicLiquid(
    name="EMI-BF4",
    conductivity=1.46,
    surface_tension=0.072,
    viscosity=0.0011,
    density=998.0,
    permittivity=80.0,
    M_cation=111.17e-3,
    M_anion=86.81e-3
)

EAN = IonicLiquid(
    name="EAN",
    conductivity=2.27,
    surface_tension=0.0481,
    viscosity=0.0398,
    density=1210.0,
    permittivity=29.0,
    M_cation=46.1e-3,
    M_anion=62.0e-3
)

BMI_TCM = IonicLiquid(
    name="BMI-TCM",
    conductivity=1.03,
    surface_tension=0.0496,
    viscosity=0.0257,
    density=1050.0,
    permittivity=14.0,
    M_cation=139.2e-3,
    M_anion=90.1e-3
)


class ElectrosprayPhysics:
    """Core physics calculations for electrospray"""
    
    def __init__(self, liquid: IonicLiquid):
        self.liquid = liquid
        self.T_ambient = 298.15  # K
        
    def characteristic_length(self, Q: float) -> float:
        """
        Gañán-Calvo characteristic length scale r_G
        
        r_G = (ρε₀Q³/γK)^(1/6)
        """
        return (self.liquid.density * EPSILON_0 * Q**3 / 
                (self.liquid.surface_tension * self.liquid.conductivity))**(1/6)
    
    def dimensionless_flow_rate(self, Q: float) -> float:
        """
        Dimensionless flow rate Π
        
        Π = ρKQ/(γε₀)
        """
        return (self.liquid.density * self.liquid.conductivity * Q / 
                (self.liquid.surface_tension * EPSILON_0))
    
    def electrohydrodynamic_reynolds(self) -> float:
        """
        Electrohydrodynamic Reynolds number Re_K
        
        Re_K = (γ²ρε₀/μ³K)^(1/3)
        """
        return ((self.liquid.surface_tension**2 * self.liquid.density * EPSILON_0) / 
                (self.liquid.viscosity**3 * self.liquid.conductivity))**(1/3)
    
    def effective_pressure_drop(self) -> float:
        """
        Effective pressure drop for electrospray (Gañán-Calvo)
        
        ΔP = k_p(σ²K²/ε₀ρ²)^(1/3)
        
        where k_p ≈ 1.3 from experiments
        """
        k_p = 1.3
        return k_p * ((self.liquid.surface_tension**2 * self.liquid.conductivity**2) / 
                      (EPSILON_0 * self.liquid.density**2))**(1/3)
    
    def jet_radius(self, Q: float, Delta_P: float = None) -> float:
        """
        Jet radius from energy balance
        
        R = (ρQ²/2π²ΔP)^(1/4)
        """
        if Delta_P is None:
            Delta_P = self.effective_pressure_drop()
        
        return (self.liquid.density * Q**2 / (2 * np.pi**2 * Delta_P))**(1/4)
    
    def emitted_current_isothermal(self, Q: float, V: float = None, d_tip: float = 30e-6) -> float:
        """
        Isothermal current emission with voltage dependence
        
        I = [I₀ + ψ√(γKṁ/ρ)] · [1 - exp(-(V - V_onset)/V_char)]
        
        Based on Gamero-Castaño (2008) formulation with Gañán-Calvo scaling
        
        Args:
            Q: Volumetric flow rate (m³/s)
            V: Applied voltage (V). If None, returns asymptotic current
            d_tip: Tip diameter (m)
        
        References:
            - Gamero-Castaño & Hruby (2002) PRE 65, 021402 (onset voltage)
            - Gamero-Castaño (2008) APL 92, 163502 (voltage dependence)
        """
        # Convert Q to mass flow rate
        m_dot = Q * self.liquid.density
        
        # Empirical parameters (from Paper 3, Figure 4)
        psi_dict = {
            'EMI-Im': (2.48, 81e-9),
            'EMI-BF4': (2.48, 81e-9),  # Approximate
            'BMI-TCM': (2.41, 217e-9),
            'EAN': (2.26, 340e-9)
        }
        
        psi, I_0 = psi_dict.get(self.liquid.name, (2.6, 0))
        
        # Asymptotic current (high voltage limit)
        I_asymptotic = I_0 + psi * np.sqrt(self.liquid.surface_tension * 
                                           self.liquid.conductivity * m_dot / 
                                           self.liquid.density)
        
        # If voltage not specified, return asymptotic value
        if V is None:
            return I_asymptotic
        
        # Calculate onset voltage (Taylor cone formation threshold)
        # V_onset = A · √(γd_tip/(ε₀ε_r))
        A_onset = 1.8  # Geometric factor from experiments (Gamero-Castaño 2002)
        V_onset = A_onset * np.sqrt(self.liquid.surface_tension * d_tip / 
                                    (EPSILON_0 * self.liquid.permittivity))
        
        # Voltage-dependent factor
        # Characteristic voltage for transition (from experiments)
        V_char = 300.0  # Volts
        
        if V < V_onset:
            return 0.0  # No emission below onset
        else:
            # Smooth onset with exponential saturation
            g_V = 1.0 - np.exp(-(V - V_onset) / V_char)
            return I_asymptotic * g_V
    
    def self_heating_temperature_rise(self, m_dot: float, location: str = 'crossover') -> float:
        """
        Temperature rise due to self-heating (Magnani & Gamero-Castaño)
        
        ΔT = b₁ṁ^(-b₂) + b₃
        
        Locations: 'crossover' or 'jet_500rG'
        """
        # Coefficients from Paper 3, Table 2
        coeffs = {
            'EMI-Im': {
                'crossover': (0.526, 0.400, -7.63),
                'jet_500rG': (0.0104, 0.390, -3.61)
            },
            'EAN': {
                'crossover': (0.479, 0.243, -39.2),
                'jet_500rG': (1.27, 0.223, -52.9)
            },
            'BMI-TCM': {
                'crossover': (0.213, 0.259, -28.2),
                'jet_500rG': (0.130, 0.292, -10.2)
            }
        }
        
        if self.liquid.name in coeffs:
            b1, b2, b3 = coeffs[self.liquid.name][location]
            return b1 * m_dot**(-b2) + b3
        else:
            # Return conservative estimate
            return 200.0  # Kelvin
    
    def emitted_current_with_heating(self, Q: float, V: float = None, d_tip: float = 30e-6) -> float:
        """
        Current emission accounting for self-heating and voltage dependence
        
        Uses temperature-dependent conductivity combined with voltage-dependent emission
        
        Args:
            Q: Volumetric flow rate (m³/s)
            V: Applied voltage (V)
            d_tip: Tip diameter (m)
        """
        m_dot = Q * self.liquid.density
        
        # Get temperature rise
        Delta_T = self.self_heating_temperature_rise(m_dot, 'crossover')
        T_eff = self.T_ambient + Delta_T
        
        # Update conductivity for temperature
        K_eff = self.liquid.conductivity_T(T_eff)
        
        # Calculate asymptotic current with effective conductivity
        psi = 2.5  # Average value
        I_0 = 100e-9  # Average intercept
        
        I_asymptotic = I_0 + psi * np.sqrt(self.liquid.surface_tension * K_eff * m_dot / 
                                           self.liquid.density)
        
        # If voltage not specified, return asymptotic value
        if V is None:
            return I_asymptotic
        
        # Calculate onset voltage
        A_onset = 1.8
        V_onset = A_onset * np.sqrt(self.liquid.surface_tension * d_tip / 
                                    (EPSILON_0 * self.liquid.permittivity))
        
        # Voltage-dependent factor
        V_char = 300.0  # Volts
        
        if V < V_onset:
            return 0.0
        else:
            g_V = 1.0 - np.exp(-(V - V_onset) / V_char)
            return I_asymptotic * g_V
    
    def minimum_flow_rate_surface_tension(self) -> float:
        """
        Minimum flow rate limited by surface tension (We = 1)
        
        Q_σ = (σ⁴/ΔP³ρ²)^(1/2)
        """
        Delta_P = self.effective_pressure_drop()
        return np.sqrt(self.liquid.surface_tension**4 / 
                      (Delta_P**3 * self.liquid.density**2))
    
    def minimum_flow_rate_viscous(self) -> float:
        """
        Minimum flow rate limited by viscosity
        
        Q_μ = (σ⁴/μ³ΔP)^(1/2)
        """
        Delta_P = self.effective_pressure_drop()
        return np.sqrt(self.liquid.surface_tension**4 / 
                      (self.liquid.viscosity**3 * Delta_P))
    
    def minimum_flow_rate(self) -> float:
        """Overall minimum stable flow rate"""
        Q_sigma = self.minimum_flow_rate_surface_tension()
        Q_mu = self.minimum_flow_rate_viscous()
        
        # Take the larger of the two
        return max(Q_sigma, Q_mu)
    
    def ion_solvation_energy(self, n: int, charge: int = 1) -> float:
        """
        Ion solvation energy barrier (Born model)
        
        ΔG°_s ∝ (ne)^(4/3)
        
        Args:
            n: solvation number
            charge: charge state
        """
        # Typical value ~1-2 eV for n=0
        Delta_G_0 = 1.5 * E_CHARGE  # J
        
        return Delta_G_0 * ((n + 1) * charge)**(4/3) / charge**(4/3)
    
    def ion_evaporation_field_reduction(self, E: float, charge: int = 1) -> float:
        """
        Field-induced reduction in solvation energy
        
        ΔGₑ = √((ne)³E/4πε₀)
        """
        n_eff = 0  # For bare ion
        return np.sqrt((charge * E_CHARGE)**3 * E / (4 * np.pi * EPSILON_0))
    
    def ion_emission_fraction(self, E: float, T: float, charge: int = 1) -> float:
        """
        Fraction of ions vs droplets (Boltzmann factor)
        
        P_ion ∝ exp(-ΔG/k_B T)
        """
        Delta_G_s = self.ion_solvation_energy(0, charge)
        Delta_G_e = self.ion_evaporation_field_reduction(E, charge)
        Delta_G = Delta_G_s - Delta_G_e
        
        return np.exp(-Delta_G / (K_B * T))


class SingleCapillaryEmitter:
    """Model for single capillary electrospray emitter"""
    
    def __init__(self, 
                 liquid: IonicLiquid,
                 d_tip: float = 30e-6,           # Tip diameter (m)
                 d_extractor: float = 1e-3,       # Extractor aperture (m)
                 t_extractor: float = 0.5e-3,     # Extractor thickness (m)
                 gap: float = 0.5e-3,             # Tip-extractor gap (m)
                 V_emitter: float = 2000.0,       # Emitter voltage (V)
                 m_dot: float = 1e-10):           # Mass flow rate (kg/s)
        
        self.liquid = liquid
        self.physics = ElectrosprayPhysics(liquid)
        
        # Geometry
        self.d_tip = d_tip
        self.d_extractor = d_extractor
        self.t_extractor = t_extractor
        self.gap = gap
        
        # Operating conditions
        self.V_emitter = V_emitter
        self.m_dot = m_dot
        self.Q = m_dot / liquid.density
        
        # Calculated parameters
        self.I = None
        self.R_jet = None
        self.R_droplet = None
        self.thrust = None
        self.Isp = None
        self.efficiency = None
        self.emission_mode = None
        
    def electric_field_at_tip(self) -> float:
        """
        Estimate electric field at emitter tip
        
        Simple approximation: E ≈ 2V/d (for hemisphere-plane geometry)
        More accurate: numerical solution of Laplace equation
        """
        # Simplified model
        E_simple = 2 * self.V_emitter / self.gap
        
        # Enhancement factor due to sharp tip
        # β ≈ 2/r_tip for hemisphere
        beta = 2 / (self.d_tip / 2)
        
        # Effective field
        E_eff = E_simple * (1 + 0.5 * np.log(beta))
        
        return E_eff
    
    def onset_voltage(self) -> float:
        """
        Calculate onset voltage for Taylor cone formation
        
        V_onset = A · √(γd_tip/(ε₀ε_r))
        
        Based on Taylor cone theory and experiments by Gamero-Castaño & Hruby (2002)
        
        Below this voltage, no stable cone-jet emission occurs.
        
        Returns:
            V_onset (V): Minimum voltage for stable emission
        """
        A_onset = 1.8  # Geometric factor from experiments
        
        V_onset = A_onset * np.sqrt(self.liquid.surface_tension * self.d_tip / 
                                    (EPSILON_0 * self.liquid.permittivity))
        
        return V_onset
    
    def voltage_for_flow_rate(self, m_dot: float) -> float:
        """
        Calculate required voltage for stable emission at given flow rate
        
        V = V_onset + β · ṁ^0.45
        
        Based on Gañán-Calvo (1997) scaling + Lozano (2005) empirical fits
        
        Args:
            m_dot: Mass flow rate (kg/s)
            
        Returns:
            V_required (V): Voltage needed for stable emission
            
        References:
            - Gañán-Calvo (1997) PRL 79, 217
            - Gamero-Castaño & Hruby (2002) PRE 65, 021402
            - Lozano & Martínez-Sánchez (2005) JCIS 282, 415
        """
        # Onset voltage (minimum)
        V_onset = self.onset_voltage()
        
        # Empirical scaling coefficient (depends on liquid & geometry)
        # Calibrated to match Gamero-Castaño (2002) data
        # At ṁ = 5 ng/s = 5e-12 kg/s → V ≈ 1900 V
        # At ṁ = 10 ng/s = 10e-12 kg/s → V ≈ 2300 V
        # Units: V·s^0.45/kg^0.45
        beta_dict = {
            'EMI-Im': 2.15e8,      # Calibrated to experiments
            'EMI-BF4': 2.0e8,
            'EAN': 1.8e8,
            'BMI-TCM': 2.3e8
        }
        beta = beta_dict.get(self.liquid.name, 2.15e8)
        
        # Required voltage
        V_required = V_onset + beta * m_dot**0.45
        
        return V_required
    
    def flow_rate_for_voltage(self, V: float) -> float:
        """
        Calculate achievable flow rate at given voltage
        
        Inverse of voltage_for_flow_rate:
        ṁ = [(V - V_onset) / β]^(1/0.45)
        
        Args:
            V: Applied voltage (V)
            
        Returns:
            m_dot (kg/s): Achievable mass flow rate
        """
        V_onset = self.onset_voltage()
        
        if V <= V_onset:
            return 0.0  # No emission below onset
        
        beta_dict = {
            'EMI-Im': 2.15e8,
            'EMI-BF4': 2.0e8,
            'EAN': 1.8e8,
            'BMI-TCM': 2.3e8
        }
        beta = beta_dict.get(self.liquid.name, 2.15e8)
        
        # Solve for flow rate
        m_dot = ((V - V_onset) / beta)**(1.0 / 0.45)
        
        return m_dot
    
    def generate_operating_envelope(self, m_dot_range: np.ndarray) -> Dict:
        """
        Generate V-ṁ operating envelope showing stable region
        
        Returns three curves:
        - V_min: Minimum voltage (onset + 10% margin)
        - V_nominal: Nominal operating voltage (stable cone-jet)
        - V_max: Maximum voltage (before multi-jet/corona)
        
        Args:
            m_dot_range: Array of mass flow rates (kg/s)
            
        Returns:
            Dict with 'flow_rates', 'V_min', 'V_nominal', 'V_max'
        """
        V_onset = self.onset_voltage()
        
        V_min = []
        V_nominal = []
        V_max = []
        
        beta_dict = {
            'EMI-Im': 2.15e8,
            'EMI-BF4': 2.0e8,
            'EAN': 1.8e8,
            'BMI-TCM': 2.3e8
        }
        beta = beta_dict.get(self.liquid.name, 2.15e8)
        
        for m_dot in m_dot_range:
            # Minimum: onset + 10% margin
            V_min.append(V_onset * 1.1)
            
            # Nominal: Gañán-Calvo scaling
            V_nom = V_onset + beta * m_dot**0.45
            V_nominal.append(V_nom)
            
            # Maximum: 1.5× nominal (before instability)
            V_max.append(V_nom * 1.5)
        
        return {
            'flow_rates': m_dot_range,
            'V_min': np.array(V_min),
            'V_nominal': np.array(V_nominal),
            'V_max': np.array(V_max),
            'V_onset': V_onset
        }
    
    def calculate_emission(self, use_heating: bool = True):
        """Calculate all emission parameters with voltage-dependent current"""
        
        # Check if flow rate is above minimum
        Q_min = self.physics.minimum_flow_rate()
        if self.Q < Q_min:
            print(f"Warning: Flow rate {self.Q:.2e} m³/s below minimum {Q_min:.2e} m³/s")
            print("Emission may be unstable!")
        
        # Calculate current with voltage dependence
        if use_heating:
            self.I = self.physics.emitted_current_with_heating(self.Q, self.V_emitter, self.d_tip)
        else:
            self.I = self.physics.emitted_current_isothermal(self.Q, self.V_emitter, self.d_tip)
        
        # Calculate jet radius
        self.R_jet = self.physics.jet_radius(self.Q)
        
        # Calculate droplet radius (≈ 1.89 × R_jet from Rayleigh breakup)
        self.R_droplet = 1.89 * self.R_jet
        
        # Determine emission mode
        self._determine_emission_mode()
        
        # Calculate performance
        self._calculate_performance()
        
        return self
    
    def _determine_emission_mode(self):
        """
        Determine if emission is ion-dominated or droplet-dominated
        
        Based on Paper 3: pure ion mode at ṁ < 4×10⁻¹² kg/s
        """
        if self.m_dot < 4e-12:
            self.emission_mode = "Pure Ion"
            self.mode_fraction_ions = 1.0
        elif self.m_dot < 1e-10:
            self.emission_mode = "Mixed (Ion-Droplet)"
            # Empirical transition
            self.mode_fraction_ions = np.exp(-(self.m_dot - 4e-12) / 2e-11)
        else:
            self.emission_mode = "Droplet-Dominated"
            self.mode_fraction_ions = 0.15  # ~15% ion current typical
    
    def _calculate_performance(self):
        """Calculate thrust, Isp, and efficiency"""
        
        # Average charge-to-mass ratio
        if self.emission_mode == "Pure Ion":
            # Typical ion cluster: A⁺[A⁺B⁻]₂ (trimer)
            M_cluster = (self.liquid.M_cation + 
                        2 * (self.liquid.M_cation + self.liquid.M_anion))
            q_m_avg = E_CHARGE / M_cluster  # C/kg
        else:
            # Droplet-dominated
            q_m_avg = self.I / self.m_dot  # C/kg
        
        # Exhaust velocity
        c_exhaust = np.sqrt(2 * self.V_emitter * q_m_avg)
        
        # Thrust
        self.thrust = self.m_dot * c_exhaust
        
        # Specific impulse
        self.Isp = c_exhaust / G_0
        
        # Efficiency
        power_in = self.I * self.V_emitter
        power_beam = 0.5 * self.m_dot * c_exhaust**2
        
        if power_in > 0:
            self.efficiency = power_beam / power_in
        else:
            self.efficiency = 0
        
        # Limit efficiency to physical range
        self.efficiency = min(self.efficiency, 0.95)
    
    def print_results(self):
        """Print comprehensive results"""
        print("\n" + "="*70)
        print(f"SINGLE CAPILLARY ELECTROSPRAY SIMULATION")
        print("="*70)
        
        print(f"\nPropellant: {self.liquid.name}")
        print(f"  Conductivity: {self.liquid.conductivity:.2f} S/m")
        print(f"  Surface Tension: {self.liquid.surface_tension*1e3:.2f} mN/m")
        print(f"  Viscosity: {self.liquid.viscosity*1e3:.2f} mPa·s")
        print(f"  Density: {self.liquid.density:.0f} kg/m³")
        
        print(f"\nGeometry:")
        print(f"  Tip Diameter: {self.d_tip*1e6:.1f} μm")
        print(f"  Extractor Gap: {self.gap*1e3:.2f} mm")
        print(f"  Extractor Aperture: {self.d_extractor*1e3:.2f} mm")
        
        print(f"\nOperating Conditions:")
        print(f"  Emitter Voltage: {self.V_emitter:.0f} V")
        print(f"  Mass Flow Rate: {self.m_dot:.2e} kg/s")
        print(f"  Volumetric Flow Rate: {self.Q:.2e} m³/s")
        
        print(f"\nCharacteristic Scales:")
        r_G = self.physics.characteristic_length(self.Q)
        print(f"  Characteristic Length (r_G): {r_G*1e9:.2f} nm")
        print(f"  Dimensionless Flow (Π): {self.physics.dimensionless_flow_rate(self.Q):.2e}")
        print(f"  EHD Reynolds (Re_K): {self.physics.electrohydrodynamic_reynolds():.2f}")
        
        print(f"\nMinimum Flow Rates:")
        Q_min = self.physics.minimum_flow_rate()
        print(f"  Overall Minimum: {Q_min:.2e} m³/s ({Q_min*self.liquid.density:.2e} kg/s)")
        print(f"  Current/Minimum Ratio: {self.Q/Q_min:.2f}")
        
        print(f"\nEmission Characteristics:")
        print(f"  Jet Radius: {self.R_jet*1e9:.2f} nm")
        print(f"  Droplet Radius: {self.R_droplet*1e6:.2f} μm")
        print(f"  Emitted Current: {self.I*1e9:.2f} nA")
        print(f"  Emission Mode: {self.emission_mode}")
        print(f"  Ion Fraction: {self.mode_fraction_ions:.1%}")
        
        print(f"\nPerformance:")
        print(f"  Thrust: {self.thrust*1e6:.2f} μN")
        print(f"  Specific Impulse: {self.Isp:.0f} s")
        print(f"  Efficiency: {self.efficiency:.1%}")
        print(f"  Power: {self.I*self.V_emitter*1e6:.2f} μW")
        
        print("\n" + "="*70)
    
    def sweep_voltage(self, V_range: np.ndarray) -> Dict:
        """Sweep voltage and return performance curves"""
        results = {
            'voltage': V_range,
            'current': [],
            'thrust': [],
            'isp': [],
            'efficiency': [],
            'mode': []
        }
        
        original_V = self.V_emitter
        
        for V in V_range:
            self.V_emitter = V
            self.calculate_emission()
            
            results['current'].append(self.I)
            results['thrust'].append(self.thrust)
            results['isp'].append(self.Isp)
            results['efficiency'].append(self.efficiency)
            results['mode'].append(self.mode_fraction_ions)
        
        self.V_emitter = original_V
        
        # Convert to arrays
        for key in ['current', 'thrust', 'isp', 'efficiency', 'mode']:
            results[key] = np.array(results[key])
        
        return results
    
    def sweep_flow_rate(self, m_dot_range: np.ndarray) -> Dict:
        """Sweep mass flow rate and return performance curves"""
        results = {
            'm_dot': m_dot_range,
            'current': [],
            'thrust': [],
            'isp': [],
            'efficiency': [],
            'mode': [],
            'jet_radius': []
        }
        
        original_m_dot = self.m_dot
        
        for m_dot in m_dot_range:
            self.m_dot = m_dot
            self.Q = m_dot / self.liquid.density
            self.calculate_emission()
            
            results['current'].append(self.I)
            results['thrust'].append(self.thrust)
            results['isp'].append(self.Isp)
            results['efficiency'].append(self.efficiency)
            results['mode'].append(self.mode_fraction_ions)
            results['jet_radius'].append(self.R_jet)
        
        self.m_dot = original_m_dot
        self.Q = original_m_dot / self.liquid.density
        
        # Convert to arrays
        for key in ['current', 'thrust', 'isp', 'efficiency', 'mode', 'jet_radius']:
            results[key] = np.array(results[key])
        
        return results


class MultiCapillaryArray:
    """Model for multi-capillary electrospray array"""
    
    def __init__(self,
                 liquid: IonicLiquid,
                 n_emitters: int = 100,
                 d_tip: float = 30e-6,
                 array_spacing: float = 1e-3,
                 d_extractor: float = 10e-3,
                 t_extractor: float = 0.5e-3,
                 gap: float = 0.5e-3,
                 V_emitter: float = 2000.0,
                 m_dot_total: float = 1e-8):
        
        self.liquid = liquid
        self.n_emitters = n_emitters
        self.array_spacing = array_spacing
        
        # Mass flow rate per emitter
        m_dot_per_emitter = m_dot_total / n_emitters
        
        # Create individual emitters
        self.emitters = [
            SingleCapillaryEmitter(
                liquid=liquid,
                d_tip=d_tip,
                d_extractor=d_extractor,
                t_extractor=t_extractor,
                gap=gap,
                V_emitter=V_emitter,
                m_dot=m_dot_per_emitter
            )
            for _ in range(n_emitters)
        ]
        
        # Array-level parameters
        self.total_current = None
        self.total_thrust = None
        self.array_Isp = None
        self.array_efficiency = None
        self.uniformity = None
    
    def calculate_array_emission(self, 
                                 uniformity: float = 0.95,
                                 use_heating: bool = True):
        """
        Calculate emission from entire array
        
        Args:
            uniformity: Fraction of emitters operating identically (0-1)
                       Accounts for manufacturing variations
        """
        self.uniformity = uniformity
        
        # Calculate for representative emitter
        self.emitters[0].calculate_emission(use_heating=use_heating)
        
        # Account for non-uniformity
        n_uniform = int(self.n_emitters * uniformity)
        n_varied = self.n_emitters - n_uniform
        
        # Uniform emitters
        I_uniform = self.emitters[0].I * n_uniform
        T_uniform = self.emitters[0].thrust * n_uniform
        
        # Varied emitters (assume ±20% variation in performance)
        if n_varied > 0:
            variation = np.random.normal(1.0, 0.2, n_varied)
            variation = np.clip(variation, 0.5, 1.5)
            
            I_varied = self.emitters[0].I * np.sum(variation)
            T_varied = self.emitters[0].thrust * np.sum(variation)
        else:
            I_varied = 0
            T_varied = 0
        
        # Total performance
        self.total_current = I_uniform + I_varied
        self.total_thrust = T_uniform + T_varied
        
        # Array-level metrics
        m_dot_total = self.emitters[0].m_dot * self.n_emitters
        self.array_Isp = self.total_thrust / (m_dot_total * G_0)
        
        power_in = self.total_current * self.emitters[0].V_emitter
        power_beam = 0.5 * m_dot_total * (self.total_thrust / m_dot_total)**2
        self.array_efficiency = power_beam / power_in if power_in > 0 else 0
        
        return self
    
    def print_array_results(self):
        """Print array-level results"""
        print("\n" + "="*70)
        print(f"MULTI-CAPILLARY ARRAY SIMULATION")
        print("="*70)
        
        print(f"\nArray Configuration:")
        print(f"  Number of Emitters: {self.n_emitters}")
        print(f"  Emitter Spacing: {self.array_spacing*1e3:.2f} mm")
        print(f"  Array Uniformity: {self.uniformity:.1%}")
        
        print(f"\nSingle Emitter Performance:")
        print(f"  Current: {self.emitters[0].I*1e9:.2f} nA")
        print(f"  Thrust: {self.emitters[0].thrust*1e6:.2f} μN")
        print(f"  Mass Flow: {self.emitters[0].m_dot:.2e} kg/s")
        
        print(f"\nTotal Array Performance:")
        print(f"  Total Current: {self.total_current*1e6:.2f} μA")
        print(f"  Total Thrust: {self.total_thrust*1e3:.2f} mN")
        print(f"  Specific Impulse: {self.array_Isp:.0f} s")
        print(f"  Efficiency: {self.array_efficiency:.1%}")
        print(f"  Total Power: {self.total_current*self.emitters[0].V_emitter:.2f} W")
        
        # Calculate thrust density
        array_area = (self.array_spacing * np.sqrt(self.n_emitters))**2
        thrust_density = self.total_thrust / array_area
        print(f"  Thrust Density: {thrust_density:.2f} N/m²")
        
        print("\n" + "="*70)
    
    def optimize_performance(self, 
                           target: str = 'thrust',
                           constraint_power: float = None,
                           constraint_isp_min: float = None):
        """
        Optimize array performance
        
        Args:
            target: 'thrust', 'efficiency', or 'isp'
            constraint_power: Maximum power (W)
            constraint_isp_min: Minimum Isp (s)
        """
        
        def objective(params):
            """Objective function to minimize (negative of target)"""
            V, m_dot_total = params
            
            # Update array
            m_dot_per = m_dot_total / self.n_emitters
            for emitter in self.emitters:
                emitter.V_emitter = V
                emitter.m_dot = m_dot_per
                emitter.Q = m_dot_per / self.liquid.density
            
            # Calculate performance
            self.calculate_array_emission()
            
            # Check constraints
            if constraint_power is not None:
                power = self.total_current * V
                if power > constraint_power:
                    return 1e10  # Penalty
            
            if constraint_isp_min is not None:
                if self.array_Isp < constraint_isp_min:
                    return 1e10  # Penalty
            
            # Return negative of target (for minimization)
            if target == 'thrust':
                return -self.total_thrust
            elif target == 'efficiency':
                return -self.array_efficiency
            elif target == 'isp':
                return -self.array_Isp
            else:
                return 0
        
        # Initial guess
        V_0 = 2000.0
        m_dot_0 = 1e-8
        
        # Bounds
        bounds = [(500, 10000), (1e-10, 1e-6)]
        
        # Optimize
        result = minimize(objective, [V_0, m_dot_0], 
                         bounds=bounds, method='L-BFGS-B')
        
        if result.success:
            V_opt, m_dot_opt = result.x
            print(f"\nOptimization successful!")
            print(f"  Optimal Voltage: {V_opt:.0f} V")
            print(f"  Optimal Total Mass Flow: {m_dot_opt:.2e} kg/s")
            
            # Set to optimal and recalculate
            for emitter in self.emitters:
                emitter.V_emitter = V_opt
                emitter.m_dot = m_dot_opt / self.n_emitters
                emitter.Q = emitter.m_dot / self.liquid.density
            
            self.calculate_array_emission()
            self.print_array_results()
        else:
            print(f"\nOptimization failed: {result.message}")
        
        return result


def demonstrate_single_emitter():
    """Demonstrate single capillary emitter simulation"""
    
    print("\n" + "#"*70)
    print("# SINGLE CAPILLARY EMITTER DEMONSTRATION")
    print("#"*70)
    
    # Create emitter with EMI-Im
    emitter = SingleCapillaryEmitter(
        liquid=EMI_IM,
        d_tip=20e-6,           # 20 μm tip
        gap=0.5e-3,            # 0.5 mm gap
        V_emitter=2000,        # 2 kV
        m_dot=5e-11            # 50 pg/s
    )
    
    # Calculate and print results
    emitter.calculate_emission(use_heating=True)
    emitter.print_results()
    
    # Voltage sweep
    print("\nPerforming voltage sweep...")
    V_range = np.linspace(1000, 5000, 50)
    results_V = emitter.sweep_voltage(V_range)
    
    # Flow rate sweep
    print("Performing flow rate sweep...")
    m_dot_range = np.logspace(-12, -9, 50)
    results_Q = emitter.sweep_flow_rate(m_dot_range)
    
    return emitter, results_V, results_Q


def demonstrate_multi_capillary():
    """Demonstrate multi-capillary array simulation"""
    
    print("\n" + "#"*70)
    print("# MULTI-CAPILLARY ARRAY DEMONSTRATION")
    print("#"*70)
    
    # Create 100-emitter array
    array = MultiCapillaryArray(
        liquid=EMI_IM,
        n_emitters=100,
        d_tip=20e-6,
        array_spacing=1e-3,
        V_emitter=2000,
        m_dot_total=5e-9       # 5 ng/s total
    )
    
    # Calculate and print results
    array.calculate_array_emission(uniformity=0.95)
    array.print_array_results()
    
    # Optimize for maximum thrust
    print("\nOptimizing for maximum thrust...")
    array.optimize_performance(
        target='thrust',
        constraint_power=10.0,      # 10 W max
        constraint_isp_min=1000     # Min 1000 s Isp
    )
    
    return array


if __name__ == "__main__":
    # Run demonstrations
    emitter, results_V, results_Q = demonstrate_single_emitter()
    array = demonstrate_multi_capillary()
    
    print("\n" + "#"*70)
    print("# Simulation complete! Results stored in variables:")
    print("#   emitter: Single capillary object")
    print("#   array: Multi-capillary array object")
    print("#   results_V: Voltage sweep data")
    print("#   results_Q: Flow rate sweep data")
    print("#"*70)