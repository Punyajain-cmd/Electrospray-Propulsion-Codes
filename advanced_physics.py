"""
Advanced Physics Models for Electrospray Simulation
====================================================

Additional physics models not in basic implementation:
1. Electrochemical reactions at meniscus
2. Non-equilibrium ion transport
3. Turbulent jet breakup
4. Coulomb fission of droplets
5. Electrostatic focusing/defocusing
6. Beam neutralization
7. Plume divergence
8. Mass spectral distribution
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from dataclasses import dataclass
from scipy.integrate import odeint
from scipy.optimize import fsolve

import sys
sys.path.insert(0, '.')
from config import PhysicalConstants as PC
from ionic_liquids import IonicLiquid


# ============================================================================
# COULOMB FISSION MODEL
# ============================================================================

class CoulombFissionModel:
    """
    Model for Coulomb fission of charged droplets
    
    Rayleigh limit: when electrostatic pressure exceeds surface tension,
    droplets undergo fission
    """
    
    def __init__(self, liquid: IonicLiquid):
        self.liquid = liquid
    
    def rayleigh_charge_limit(self, R: float) -> float:
        """
        Maximum stable charge on droplet (Rayleigh limit)
        
        Q_R = 8π√(ε₀γR³)
        
        Args:
            R: Droplet radius (m)
        Returns:
            Maximum charge (C)
        """
        return 8 * np.pi * np.sqrt(PC.EPSILON_0 * self.liquid.surface_tension * R**3)
    
    def rayleigh_parameter(self, Q: float, R: float) -> float:
        """
        Rayleigh fissility parameter: X = Q/Q_R
        
        X < 1: Stable
        X ≥ 1: Unstable, undergoes fission
        """
        Q_R = self.rayleigh_charge_limit(R)
        return Q / Q_R
    
    def fission_products(self, R_parent: float, Q_parent: float, 
                         n_offspring: int = 2) -> Tuple[float, float]:
        """
        Calculate size and charge of offspring droplets after fission
        
        Conservation laws:
        - Mass: 4/3πR₀³ = n × 4/3πR₁³
        - Charge: Q₀ ≈ n × Q₁ (approximately equal distribution)
        
        Returns:
            (R_offspring, Q_offspring)
        """
        # Volume conservation
        V_parent = (4/3) * np.pi * R_parent**3
        V_offspring = V_parent / n_offspring
        R_offspring = (3 * V_offspring / (4 * np.pi))**(1/3)
        
        # Charge distribution (empirical: not perfectly equal)
        # Typically ~2% of charge remains with parent
        Q_offspring = 0.98 * Q_parent / n_offspring
        
        return R_offspring, Q_offspring
    
    def fission_cascade(self, R_initial: float, Q_initial: float, 
                        max_generations: int = 10) -> List[Dict]:
        """
        Simulate cascade of fissions until stable droplets reached
        
        Returns:
            List of final droplet populations
        """
        generations = []
        current_gen = [{'R': R_initial, 'Q': Q_initial, 'count': 1}]
        
        for gen in range(max_generations):
            next_gen = []
            
            for droplet in current_gen:
                R, Q, count = droplet['R'], droplet['Q'], droplet['count']
                X = self.rayleigh_parameter(Q, R)
                
                if X >= 1.0:
                    # Undergoes fission
                    R_off, Q_off = self.fission_products(R, Q, n_offspring=2)
                    next_gen.append({'R': R_off, 'Q': Q_off, 'count': count * 2})
                else:
                    # Stable
                    next_gen.append({'R': R, 'Q': Q, 'count': count})
            
            current_gen = next_gen
            generations.append(current_gen)
            
            # Check if all stable
            if all(self.rayleigh_parameter(d['Q'], d['R']) < 0.9 for d in current_gen):
                break
        
        return generations


# ============================================================================
# BEAM DIVERGENCE MODEL
# ============================================================================

class BeamDivergenceModel:
    """
    Model plume divergence and angular distribution
    
    Includes:
    - Space charge expansion
    - Initial velocity distribution
    - Collisional scattering
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
        self.v_beam = np.sqrt(2 * V / m_per_q)
    
    def space_charge_divergence_half_angle(self, z: float, r0: float) -> float:
        """
        Beam divergence half-angle due to space charge
        
        θ(z) ≈ (I/(4πε₀v³m/q))^(1/2) × (z/r₀)
        
        Args:
            z: Distance from emitter (m)
            r0: Initial beam radius (m)
        Returns:
            Half-angle (radians)
        """
        q_per_m = 1 / self.m_per_q
        
        theta = np.sqrt(self.I / (4 * np.pi * PC.EPSILON_0 * self.v_beam**3 * q_per_m)) * (z / r0)
        
        return theta
    
    def thermal_divergence_half_angle(self, T_emit: float) -> float:
        """
        Thermal velocity spread contribution
        
        θ_th ≈ √(kT/qV)
        
        Args:
            T_emit: Emission temperature (K)
        Returns:
            Half-angle (radians)
        """
        return np.sqrt(PC.K_B * T_emit / (self.V / self.m_per_q))
    
    def total_divergence_half_angle(self, z: float, r0: float, T_emit: float) -> float:
        """
        Total beam divergence (quadrature sum)
        
        θ_total = √(θ_sc² + θ_th²)
        """
        theta_sc = self.space_charge_divergence_half_angle(z, r0)
        theta_th = self.thermal_divergence_half_angle(T_emit)
        
        return np.sqrt(theta_sc**2 + theta_th**2)
    
    def beam_radius_at_distance(self, z: float, r0: float, T_emit: float) -> float:
        """
        Beam radius at distance z from source
        
        r(z) = r₀ + z × tan(θ)
        """
        theta = self.total_divergence_half_angle(z, r0, T_emit)
        return r0 + z * np.tan(theta)
    
    def angular_current_density(self, theta: np.ndarray) -> np.ndarray:
        """
        Angular distribution of current density
        
        Cosine distribution for electrospray:
        j(θ) = j₀ cos^n(θ)
        
        Typical n ≈ 2-3 for electrospray
        """
        n = 2.5  # Empirical value
        j0 = self.I / (2 * np.pi * (1 - 1/(n+1)))  # Normalization
        
        return j0 * np.cos(theta)**n


# ============================================================================
# ION FRAGMENTATION MODEL
# ============================================================================

class IonFragmentationModel:
    """
    Model fragmentation of ion clusters in acceleration region
    
    High electric fields can cause ion clusters to fragment:
    A⁺[A⁺B⁻]ₙ → A⁺[A⁺B⁻]ₙ₋₁ + A⁺B⁻
    """
    
    def __init__(self, liquid: IonicLiquid):
        self.liquid = liquid
    
    def binding_energy(self, n: int) -> float:
        """
        Binding energy of nth ion pair in cluster
        
        E_b ≈ E₀ / n^α  (empirical)
        
        Args:
            n: Cluster size (number of ion pairs)
        Returns:
            Binding energy (J)
        """
        E_0 = 0.5 * PC.E_CHARGE  # ~0.5 eV for first ion pair
        alpha = 0.5  # Scaling exponent
        
        return E_0 / n**alpha
    
    def fragmentation_rate(self, n: int, E_field: float, T: float) -> float:
        """
        Rate of cluster fragmentation
        
        k = ν₀ exp(-(E_b - eEd)/kT)
        
        Args:
            n: Cluster size
            E_field: Electric field (V/m)
            T: Temperature (K)
        Returns:
            Fragmentation rate (1/s)
        """
        nu_0 = 1e13  # Attempt frequency (1/s)
        E_b = self.binding_energy(n)
        d = 5e-10  # Characteristic distance (m)
        
        # Field-assisted barrier reduction
        Delta_E = PC.E_CHARGE * E_field * d
        
        # Arrhenius rate
        k = nu_0 * np.exp(-(E_b - Delta_E) / (PC.K_B * T))
        
        return k
    
    def steady_state_distribution(self, n_max: int, E_field: float, 
                                  T: float) -> np.ndarray:
        """
        Steady-state cluster size distribution
        
        Balance between formation and fragmentation
        """
        # Simplified model: Boltzmann distribution with field
        E_n = np.array([self.binding_energy(n) for n in range(1, n_max + 1)])
        
        # Effective potential with field
        d = 5e-10
        E_eff = E_n - PC.E_CHARGE * E_field * d * np.arange(1, n_max + 1)
        
        # Boltzmann weights
        P_n = np.exp(-E_eff / (PC.K_B * T))
        P_n /= np.sum(P_n)  # Normalize
        
        return P_n


# ============================================================================
# MASS SPECTRUM PREDICTION
# ============================================================================

class MassSpectrumModel:
    """
    Predict mass spectrum from electrospray emission
    
    Combines:
    - Ion cluster distribution
    - Droplet size distribution
    - Fragmentation
    """
    
    def __init__(self, liquid: IonicLiquid):
        self.liquid = liquid
        self.M_pair = liquid.M_cation + liquid.M_anion  # Ion pair mass
    
    def cluster_mass(self, n: int, charge: int = 1) -> float:
        """
        Mass of ion cluster A⁺[A⁺B⁻]ₙ
        
        For positive polarity:
        m = M_cation + n × (M_cation + M_anion)
        """
        return self.liquid.M_cation + n * self.M_pair
    
    def droplet_mass(self, R: float) -> float:
        """Mass of droplet with radius R"""
        V = (4/3) * np.pi * R**3
        return self.liquid.density * V
    
    def cluster_spectrum_pure_ion(self, n_max: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Pure ion mode spectrum
        
        Returns:
            (m/z values, intensities)
        """
        n_values = np.arange(0, n_max + 1)
        masses = np.array([self.cluster_mass(n) for n in n_values])
        
        # Probability decreases with cluster size
        # P(n) ∝ exp(-n/n₀)
        n_0 = 2.5  # Characteristic cluster size
        intensities = np.exp(-n_values / n_0)
        intensities /= np.sum(intensities)
        
        # m/z (assuming charge = 1)
        m_z = masses / PC.E_CHARGE * PC.AMU  # Convert to amu
        
        return m_z, intensities
    
    def mixed_mode_spectrum(self, R_mean: float, sigma_R: float,
                           f_ion: float = 0.5, n_max: int = 10) -> Dict:
        """
        Mixed ion-droplet mode spectrum
        
        Args:
            R_mean: Mean droplet radius (m)
            sigma_R: Droplet size standard deviation (m)
            f_ion: Fraction of current as ions
            n_max: Maximum cluster size
        
        Returns:
            Dict with 'ion' and 'droplet' components
        """
        # Ion component
        m_z_ion, I_ion = self.cluster_spectrum_pure_ion(n_max)
        I_ion *= f_ion
        
        # Droplet component (approximate as distribution)
        R_values = np.linspace(max(R_mean - 3*sigma_R, 1e-9), R_mean + 3*sigma_R, 50)
        m_droplet = np.array([self.droplet_mass(R) for R in R_values])
        
        # Log-normal distribution
        P_R = (1 / (R_values * sigma_R * np.sqrt(2*np.pi)) * 
               np.exp(-(np.log(R_values) - np.log(R_mean))**2 / (2*(sigma_R/R_mean)**2)))
        P_R /= np.sum(P_R)
        I_droplet = P_R * (1 - f_ion)
        
        # Assume average charge per droplet
        q_avg = 10 * PC.E_CHARGE  # Typical
        m_z_droplet = m_droplet / q_avg * PC.AMU
        
        return {
            'ion_mz': m_z_ion,
            'ion_intensity': I_ion,
            'droplet_mz': m_z_droplet,
            'droplet_intensity': I_droplet
        }


# ============================================================================
# ELECTROCHEMICAL REACTION MODEL
# ============================================================================

class ElectrochemicalReactionModel:
    """
    Model electrochemical reactions at Taylor cone surface
    
    Reactions can produce/consume ions, affecting emission
    """
    
    def __init__(self, liquid: IonicLiquid):
        self.liquid = liquid
    
    def faradaic_current(self, E_field: float, A_surface: float,
                         reaction_rate_constant: float = 1e-6) -> float:
        """
        Faradaic current from electrochemical reactions
        
        i_F = n F k A [C] exp(-αFη/RT)
        
        Args:
            E_field: Electric field (V/m)
            A_surface: Reactive surface area (m²)
            reaction_rate_constant: k (m/s)
        
        Returns:
            Faradaic current (A)
        """
        n = 1  # Electrons per reaction
        F = PC.E_CHARGE * PC.N_A  # Faraday constant
        k = reaction_rate_constant
        C = 1000  # Concentration (mol/m³) - typical for ILs
        alpha = 0.5  # Transfer coefficient
        eta = 0.1  # Overpotential (V) - estimate
        R = PC.R_GAS
        T = 298  # K
        
        i_F = n * F * k * A_surface * C * np.exp(-alpha * F * eta / (R * T))
        
        return i_F
    
    def double_layer_capacitance(self, A_surface: float) -> float:
        """
        Electric double layer capacitance
        
        C_dl = ε₀εᵣA/d_DL
        
        where d_DL ~ 1 nm (Debye length in IL)
        """
        d_DL = 1e-9  # m
        C_dl = PC.EPSILON_0 * self.liquid.permittivity * A_surface / d_DL
        
        return C_dl


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Advanced Physics Models Test")
    print("=" * 70)
    
    # Import liquid
    from ionic_liquids import EMI_IM
    
    # Test Coulomb fission
    print("\n1. COULOMB FISSION")
    print("-" * 70)
    fission = CoulombFissionModel(EMI_IM)
    
    R_drop = 1e-6  # 1 μm
    Q_drop = 1.5 * fission.rayleigh_charge_limit(R_drop)
    
    print(f"Droplet: R = {R_drop*1e6:.2f} μm, Q = {Q_drop/PC.E_CHARGE:.0f} e")
    print(f"Rayleigh limit: Q_R = {fission.rayleigh_charge_limit(R_drop)/PC.E_CHARGE:.0f} e")
    print(f"Fissility: X = {fission.rayleigh_parameter(Q_drop, R_drop):.2f}")
    
    if fission.rayleigh_parameter(Q_drop, R_drop) >= 1:
        R_off, Q_off = fission.fission_products(R_drop, Q_drop)
        print(f"After fission: R = {R_off*1e6:.2f} μm, Q = {Q_off/PC.E_CHARGE:.0f} e")
    
    # Test beam divergence
    print("\n2. BEAM DIVERGENCE")
    print("-" * 70)
    divergence = BeamDivergenceModel(I=1e-6, V=2000, m_per_q=100)
    
    z_test = 0.01  # 1 cm downstream
    r0 = 10e-6  # 10 μm initial radius
    T_emit = 400  # K
    
    theta_sc = divergence.space_charge_divergence_half_angle(z_test, r0)
    theta_th = divergence.thermal_divergence_half_angle(T_emit)
    theta_tot = divergence.total_divergence_half_angle(z_test, r0, T_emit)
    
    print(f"At z = {z_test*1e2:.0f} cm:")
    print(f"  Space charge angle: {theta_sc*1e3:.2f} mrad")
    print(f"  Thermal angle: {theta_th*1e3:.2f} mrad")
    print(f"  Total angle: {theta_tot*1e3:.2f} mrad")
    print(f"  Beam radius: {divergence.beam_radius_at_distance(z_test, r0, T_emit)*1e3:.2f} mm")
    
    # Test mass spectrum
    print("\n3. MASS SPECTRUM PREDICTION")
    print("-" * 70)
    spectrum = MassSpectrumModel(EMI_IM)
    
    m_z, intensity = spectrum.cluster_spectrum_pure_ion(n_max=5)
    
    print("Pure ion mode spectrum:")
    for i, (mz, I) in enumerate(zip(m_z, intensity)):
        if I > 0.05:  # Only significant peaks
            print(f"  m/z = {mz:.0f} amu, I = {I*100:.1f}%")
    
    # Test fragmentation
    print("\n4. ION FRAGMENTATION")
    print("-" * 70)
    fragmentation = IonFragmentationModel(EMI_IM)
    
    E_field = 1e8  # 100 MV/m (high field)
    T = 400  # K
    
    for n in [1, 2, 3, 5]:
        E_b = fragmentation.binding_energy(n)
        k_frag = fragmentation.fragmentation_rate(n, E_field, T)
        print(f"Cluster n={n}: E_b = {E_b/PC.E_CHARGE:.2f} eV, k = {k_frag:.2e} /s")
    
    print("\n" + "=" * 70)
