"""
Enhanced Electrospray Physics Models
=====================================

Comprehensive physics models including:
- Core Gañán-Calvo scaling laws
- Self-heating effects  
- Ion evaporation and cluster formation
- Polydispersity modeling
- Space charge effects
- Instability analysis
"""

import numpy as np
from typing import Tuple, Optional, Dict
from dataclasses import dataclass

from config import PhysicalConstants as PC
from ionic_liquids import IonicLiquid


# ============================================================================
# CORE PHYSICS CLASS
# ============================================================================

class ElectrosprayPhysics:
    """
    Core physics calculations for electrospray emission
    
    Based on:
    - Gañán-Calvo scaling theory
    - Magnani & Gamero-Castaño self-heating model
    - Ion evaporation models
    """
    
    def __init__(self, liquid: IonicLiquid, T_ambient: float = 298.15):
        self.liquid = liquid
        self.T_ambient = T_ambient
        
    # ========================================================================
    # CHARACTERISTIC SCALES
    # ========================================================================
    
    def characteristic_length(self, Q: float) -> float:
        """
        Gañán-Calvo characteristic length scale r_G
        
        r_G = (ρε₀Q³/γK)^(1/6)
        
        This is the natural length scale for electrospray emission
        """
        return (self.liquid.density * PC.EPSILON_0 * Q**3 / 
                (self.liquid.surface_tension * self.liquid.conductivity))**(1/6)
    
    def dimensionless_flow_rate(self, Q: float) -> float:
        """
        Dimensionless flow rate Π
        
        Π = ρKQ/(γε₀)
        
        Characterizes the flow regime
        """
        return (self.liquid.density * self.liquid.conductivity * Q / 
                (self.liquid.surface_tension * PC.EPSILON_0))
    
    def electrohydrodynamic_reynolds(self) -> float:
        """
        Electrohydrodynamic Reynolds number Re_K
        
        Re_K = (γ²ρε₀/μ³K)^(1/3)
        
        Ratio of electric to viscous forces
        """
        return ((self.liquid.surface_tension**2 * self.liquid.density * PC.EPSILON_0) / 
                (self.liquid.viscosity**3 * self.liquid.conductivity))**(1/3)
    
    def weber_number(self, Q: float, R: float) -> float:
        """
        Weber number: ratio of inertial to surface tension forces
        
        We = ρQ²/(2πR³γ)
        
        We < 1 for stable emission
        """
        return (self.liquid.density * Q**2 / 
                (2 * np.pi * R**3 * self.liquid.surface_tension))
    
    # ========================================================================
    # PRESSURE AND FLOW
    # ========================================================================
    
    def effective_pressure_drop(self) -> float:
        """
        Effective pressure drop for electrospray (Gañán-Calvo)
        
        ΔP = k_p(σ²K²/ε₀ρ²)^(1/3)
        
        where k_p ≈ 1.3 from experiments
        """
        k_p = 1.3
        return k_p * ((self.liquid.surface_tension**2 * self.liquid.conductivity**2) / 
                      (PC.EPSILON_0 * self.liquid.density**2))**(1/3)
    
    def jet_radius(self, Q: float, Delta_P: Optional[float] = None) -> float:
        """
        Jet radius from energy balance
        
        R = (ρQ²/2π²ΔP)^(1/4)
        """
        if Delta_P is None:
            Delta_P = self.effective_pressure_drop()
        
        return (self.liquid.density * Q**2 / (2 * np.pi**2 * Delta_P))**(1/4)
    
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
        return max(Q_sigma, Q_mu)
    
    # ========================================================================
    # CURRENT EMISSION
    # ========================================================================
    
    def emitted_current_isothermal(self, Q: float) -> float:
        """
        Isothermal current emission (Gañán-Calvo scaling with intercept)
        
        I = I₀ + ψ√(γKṁ/ρ)
        
        Fitting parameters from experimental data
        """
        m_dot = Q * self.liquid.density
        
        # Empirical parameters (from Paper 3, Figure 4)
        psi_dict = {
            'EMI-Im': (2.48, 81e-9),
            'EMI-BF4': (2.48, 81e-9),
            'BMI-TCM': (2.41, 217e-9),
            'EAN': (2.26, 340e-9)
        }
        
        psi, I_0 = psi_dict.get(self.liquid.name, (2.6, 100e-9))
        
        I = I_0 + psi * np.sqrt(self.liquid.surface_tension * 
                                 self.liquid.conductivity * m_dot / 
                                 self.liquid.density)
        
        return I
    
    def emitted_current_with_heating(self, Q: float) -> float:
        """
        Current emission accounting for self-heating
        
        Uses temperature-dependent conductivity
        """
        m_dot = Q * self.liquid.density
        
        # Get temperature rise
        Delta_T = self.self_heating_temperature_rise(m_dot, 'crossover')
        T_eff = self.T_ambient + Delta_T
        
        # Ensure temperature is physical
        T_eff = max(T_eff, self.T_ambient)
        T_eff = min(T_eff, 600)  # Cap at 600 K
        
        # Update conductivity for temperature
        K_eff = self.liquid.conductivity_T(T_eff)
        
        # Calculate current with effective conductivity
        psi = 2.5  # Average value
        I_0 = 100e-9  # Average intercept
        
        I = I_0 + psi * np.sqrt(self.liquid.surface_tension * K_eff * m_dot / 
                                 self.liquid.density)
        
        return I
    
    # ========================================================================
    # SELF-HEATING MODEL
    # ========================================================================
    
    def self_heating_temperature_rise(self, m_dot: float, 
                                      location: str = 'crossover') -> float:
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
            Delta_T = b1 * m_dot**(-b2) + b3
        else:
            # Conservative estimate for unknown liquids
            # Assume moderate heating
            Delta_T = 200.0 * (1e-10 / max(m_dot, 1e-15))**0.3
        
        return Delta_T
    
    # ========================================================================
    # ION EVAPORATION AND CLUSTER FORMATION
    # ========================================================================
    
    def ion_solvation_energy(self, n: int, charge: int = 1) -> float:
        """
        Ion solvation energy barrier (Born model)
        
        ΔG°_s ∝ (ne)^(4/3)
        
        Args:
            n: solvation number
            charge: charge state
        """
        Delta_G_0 = 1.5 * PC.E_CHARGE  # J (typical ~1-2 eV for n=0)
        return Delta_G_0 * ((n + 1) * charge)**(4/3) / charge**(4/3)
    
    def ion_evaporation_field_reduction(self, E: float, charge: int = 1) -> float:
        """
        Field-induced reduction in solvation energy
        
        ΔGₑ = √((ne)³E/4πε₀)
        """
        n_eff = 0  # For bare ion
        return np.sqrt((charge * PC.E_CHARGE)**3 * E / (4 * np.pi * PC.EPSILON_0))
    
    def ion_emission_fraction(self, E: float, T: float, charge: int = 1) -> float:
        """
        Fraction of ions vs droplets (Boltzmann factor)
        
        P_ion ∝ exp(-ΔG/k_B T)
        """
        Delta_G_s = self.ion_solvation_energy(0, charge)
        Delta_G_e = self.ion_evaporation_field_reduction(E, charge)
        Delta_G = Delta_G_s - Delta_G_e
        
        return np.exp(-Delta_G / (PC.K_B * T))
    
    def cluster_size_distribution(self, n_max: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Distribution of cluster sizes (solvation numbers)
        
        Returns: (n_values, probabilities)
        """
        n_values = np.arange(0, n_max + 1)
        
        # Energy landscape favors small clusters
        # P(n) ∝ exp(-E(n)/kT)
        E_n = np.array([self.ion_solvation_energy(n) for n in n_values])
        P_n = np.exp(-E_n / (PC.K_B * self.T_ambient))
        
        # Normalize
        P_n /= np.sum(P_n)
        
        return n_values, P_n


# ============================================================================
# POLYDISPERSITY MODEL
# ============================================================================

class PolydispersityModel:
    """
    Model droplet size distribution (polydispersity)
    
    Droplets formed by Rayleigh breakup have a distribution of sizes
    """
    
    def __init__(self, R_mean: float, sigma_relative: float = 0.15):
        """
        Args:
            R_mean: Mean droplet radius
            sigma_relative: Relative standard deviation (σ/R_mean)
        """
        self.R_mean = R_mean
        self.sigma = sigma_relative * R_mean
    
    def size_distribution(self, R_values: np.ndarray) -> np.ndarray:
        """
        Log-normal distribution of droplet radii
        
        P(R) = (1/Rσ√2π) exp(-(ln R - ln R₀)²/2σ²)
        """
        # Log-normal parameters
        mu = np.log(self.R_mean) - 0.5 * (self.sigma / self.R_mean)**2
        sigma_log = self.sigma / self.R_mean
        
        P_R = (1 / (R_values * sigma_log * np.sqrt(2 * np.pi)) * 
               np.exp(-(np.log(R_values) - mu)**2 / (2 * sigma_log**2)))
        
        return P_R
    
    def mass_weighted_distribution(self, R_values: np.ndarray) -> np.ndarray:
        """Mass-weighted distribution (proportional to R³)"""
        P_R = self.size_distribution(R_values)
        P_m = P_R * R_values**3
        P_m /= np.trapz(P_m, R_values)  # Normalize
        return P_m
    
    def sauter_mean_diameter(self) -> float:
        """
        Sauter mean diameter: D₃₂ = Σ(n_i D_i³) / Σ(n_i D_i²)
        
        For log-normal: D₃₂ = D₀ exp(5σ²/2)
        """
        sigma_log = self.sigma / self.R_mean
        return 2 * self.R_mean * np.exp(2.5 * sigma_log**2)


# ============================================================================
# SPACE CHARGE EFFECTS
# ============================================================================

class SpaceChargeModel:
    """
    Model space charge effects on beam expansion and current limits
    
    Space charge causes beam divergence and limits current density
    """
    
    def __init__(self, I: float, V: float, m_per_q: float):
        """
        Args:
            I: Beam current (A)
            V: Acceleration voltage (V)
            m_per_q: Mass-to-charge ratio (kg/C)
        """
        self.I = I
        self.V = V
        self.m_per_q = m_per_q
    
    def child_langmuir_limit(self, gap: float, area: float) -> float:
        """
        Child-Langmuir space charge limited current
        
        I_CL = (4ε₀/9) √(2q/m) (V^(3/2) / gap²) × area
        """
        q_per_m = 1 / self.m_per_q
        I_CL = (4 * PC.EPSILON_0 / 9) * np.sqrt(2 * q_per_m) * (
            self.V**(3/2) / gap**2) * area
        
        return I_CL
    
    def beam_perveance(self) -> float:
        """
        Beam perveance: P = I / V^(3/2)
        
        Characterizes space charge strength
        """
        return self.I / self.V**(3/2)
    
    def space_charge_expansion_angle(self, beam_radius: float, 
                                     propagation_distance: float) -> float:
        """
        Beam divergence half-angle due to space charge
        
        θ ≈ √(I/(4πε₀V)) × (z/r₀)
        
        Returns: angle in radians
        """
        q_per_m = 1 / self.m_per_q
        theta = np.sqrt(self.I / (4 * np.pi * PC.EPSILON_0 * 
                                  np.sqrt(2 * q_per_m * self.V))) * (
            propagation_distance / beam_radius)
        
        return theta


# ============================================================================
# INSTABILITY ANALYSIS
# ============================================================================

@dataclass
class InstabilityAnalysis:
    """Analyze various instability mechanisms"""
    
    liquid: IonicLiquid
    Q: float
    R_jet: float
    
    def rayleigh_instability_wavelength(self) -> float:
        """
        Most unstable wavelength for Rayleigh breakup
        
        λ_max = 9.016 × R_jet
        """
        return 9.016 * self.R_jet
    
    def rayleigh_growth_rate(self) -> float:
        """
        Maximum growth rate for Rayleigh instability
        
        ω_max = √(γ / ρR³)
        """
        return np.sqrt(self.liquid.surface_tension / 
                      (self.liquid.density * self.R_jet**3))
    
    def capillary_instability_time(self) -> float:
        """Characteristic time for capillary breakup"""
        omega = self.rayleigh_growth_rate()
        return 2 * np.pi / omega
    
    def is_stable(self) -> Tuple[bool, str]:
        """
        Check if current operating point is stable
        
        Returns: (is_stable, reason)
        """
        physics = ElectrosprayPhysics(self.liquid)
        
        # Check Weber number
        We = physics.weber_number(self.Q, self.R_jet)
        if We > 1:
            return False, f"Weber number {We:.2f} > 1: inertial instability"
        
        # Check minimum flow rate
        Q_min = physics.minimum_flow_rate()
        if self.Q < Q_min:
            return False, f"Flow rate {self.Q:.2e} < minimum {Q_min:.2e}"
        
        # Check Reynolds number
        Re_K = physics.electrohydrodynamic_reynolds()
        if Re_K < 1:
            return False, f"EHD Reynolds {Re_K:.2f} < 1: viscous instability"
        
        return True, "Operating point is stable"


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    from ionic_liquids import EMI_IM
    
    print("Enhanced Electrospray Physics Models")
    print("=" * 70)
    
    # Create physics instance
    physics = ElectrosprayPhysics(EMI_IM, T_ambient=298.15)
    
    # Test parameters
    Q = 5e-14  # m³/s
    m_dot = Q * EMI_IM.density
    
    print(f"\nTest conditions:")
    print(f"  Liquid: {EMI_IM.name}")
    print(f"  Flow rate: {Q:.2e} m³/s ({m_dot:.2e} kg/s)")
    
    # Characteristic scales
    print(f"\nCharacteristic scales:")
    r_G = physics.characteristic_length(Q)
    Pi = physics.dimensionless_flow_rate(Q)
    Re_K = physics.electrohydrodynamic_reynolds()
    print(f"  r_G = {r_G*1e9:.2f} nm")
    print(f"  Π = {Pi:.2e}")
    print(f"  Re_K = {Re_K:.2f}")
    
    # Jet properties
    print(f"\nJet properties:")
    R_jet = physics.jet_radius(Q)
    Q_min = physics.minimum_flow_rate()
    print(f"  Jet radius: {R_jet*1e9:.2f} nm")
    print(f"  Minimum Q: {Q_min:.2e} m³/s")
    print(f"  Q/Q_min: {Q/Q_min:.2f}")
    
    # Current emission
    print(f"\nCurrent emission:")
    I_iso = physics.emitted_current_isothermal(Q)
    I_heat = physics.emitted_current_with_heating(Q)
    Delta_T = physics.self_heating_temperature_rise(m_dot)
    print(f"  Isothermal: {I_iso*1e9:.2f} nA")
    print(f"  With heating: {I_heat*1e9:.2f} nA")
    print(f"  Temperature rise: {Delta_T:.1f} K")
    
    # Polydispersity
    print(f"\nPolydispersity:")
    R_droplet = 1.89 * R_jet
    poly = PolydispersityModel(R_droplet, sigma_relative=0.15)
    D32 = poly.sauter_mean_diameter()
    print(f"  Mean droplet radius: {R_droplet*1e6:.3f} μm")
    print(f"  Sauter mean diameter: {D32*1e6:.3f} μm")
    
    # Instability analysis
    print(f"\nInstability analysis:")
    instability = InstabilityAnalysis(EMI_IM, Q, R_jet)
    is_stable, reason = instability.is_stable()
    lambda_max = instability.rayleigh_instability_wavelength()
    t_break = instability.capillary_instability_time()
    print(f"  Stable: {is_stable}")
    print(f"  Reason: {reason}")
    print(f"  Rayleigh wavelength: {lambda_max*1e6:.2f} μm")
    print(f"  Breakup time: {t_break*1e6:.2f} μs")
    
    print("\n" + "=" * 70)
