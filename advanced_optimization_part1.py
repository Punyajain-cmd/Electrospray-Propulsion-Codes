"""
PART 1: Advanced Geometry Optimization with Physics-Based Coupling
"""

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.spatial import distance_matrix, Voronoi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, Polygon
import matplotlib.cm as cm
from typing import Dict, List, Tuple, Optional
import json


# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

class PhysicalConstants:
    """Universal physical constants"""
    eps0 = 8.854187817e-12  # Vacuum permittivity (F/m)
    mu0 = 1.25663706212e-6  # Vacuum permeability (H/m)
    e = 1.602176634e-19     # Elementary charge (C)
    k_B = 1.380649e-23      # Boltzmann constant (J/K)
    sigma_SB = 5.670374419e-8  # Stefan-Boltzmann (W/(m²·K⁴))
    c = 299792458           # Speed of light (m/s)
    AMU = 1.66053906660e-27  # Atomic mass unit (kg)


# ============================================================================
# EMITTER-EMITTER COUPLING MODELS
# ============================================================================

class EmitterCouplingPhysics:
    """
    Rigorous physics-based coupling between emitters
    
    Mechanisms modeled:
    1. Electric field mutual shielding/enhancement
    2. Space charge beam repulsion
    3. Plume overlap and collisional effects
    4. Thermal conduction through substrate
    5. Radiative heating between emitters
    """
    
    def __init__(self, liquid_properties: Dict, operating_conditions: Dict):
        """
        Args:
            liquid_properties: EMI-Im properties
            operating_conditions: V, I, flow rate, etc.
        """
        self.liquid = liquid_properties
        self.V = operating_conditions['voltage']
        self.I_single = operating_conditions['current_per_emitter']
        self.m_dot = operating_conditions['mass_flow_per_emitter']
        self.d_tip = operating_conditions['tip_diameter']
        
        # Derived quantities
        self.P_single = self.V * self.I_single  # Power per emitter (W)
        self.thrust_single = operating_conditions.get('thrust_per_emitter', 1e-6)  # N
        
        # Beam parameters
        self.m_ion = operating_conditions.get('ion_mass', 300) * PhysicalConstants.AMU
        self.v_ion = np.sqrt(2 * self.V * PhysicalConstants.e / self.m_ion)
        self.divergence_angle = operating_conditions.get('divergence_angle', 0.05)  # rad
        
        print(f"Emitter Coupling Physics initialized:")
        print(f"  Operating: {self.V}V, {self.I_single*1e9:.1f}nA, {self.P_single*1e6:.1f}μW")
        print(f"  Ion velocity: {self.v_ion:.0f} m/s")
        print(f"  Beam divergence: {np.degrees(self.divergence_angle):.2f}°")
    
    def electric_field_coupling(self, r: float) -> float:
        """
        Electric field mutual coupling factor
        
        Args:
            r: Distance between emitters (m)
        
        Returns:
            coupling_factor: Multiplicative factor on emission (0-1+)
        
        Physics:
        - Close range (r < 10×d_tip): Shielding (reduces field)
        - Medium range (10-50×d_tip): Enhancement (field concentration)
        - Far range (r > 50×d_tip): Independent
        """
        r_norm = r / self.d_tip
        
        if r_norm < 10:
            # Strong shielding - emitters screen each other's fields
            # Exponential decay of effectiveness
            factor = 0.6 + 0.4 * (r_norm / 10)
            
        elif r_norm < 50:
            # Moderate range - slight field enhancement
            # Like parallel plate effect
            enhancement = 0.08 * np.exp(-(r_norm - 10) / 15)
            factor = 1.0 + enhancement
            
        else:
            # Far field - independent
            factor = 1.0
        
        return factor
    
    def space_charge_repulsion(self, r: float, r_vec: np.ndarray) -> np.ndarray:
        """
        Space charge repulsion between beams
        
        Args:
            r: Distance between emitters (m)
            r_vec: Vector from emitter j to emitter i (m)
        
        Returns:
            E_sc: Space charge electric field at emitter i due to beam j (V/m)
        
        Physics:
        E_sc ~ I/(2πε₀v·r) for cylindrical beam (paraxial approximation)
        """
        if r < 1e-10:
            return np.zeros(3)
        
        # Beam is cylindrical with space charge
        # Field at distance r from axis
        E_magnitude = self.I_single / (2 * np.pi * PhysicalConstants.eps0 * self.v_ion * r)
        
        # Direction: radial from beam axis
        # Check if beams actually overlap (based on divergence)
        z_distance = abs(r_vec[2]) if len(r_vec) > 2 else 0
        transverse_distance = np.sqrt(r_vec[0]**2 + r_vec[1]**2)
        
        # Beam radius at distance z
        r_beam = z_distance * np.tan(self.divergence_angle)
        
        if transverse_distance > 2 * r_beam:
            # Beams don't overlap
            return np.zeros(3)
        
        # Overlapping beams - repulsion
        E_vec = E_magnitude * (r_vec / r)
        
        return E_vec
    
    def plume_overlap_factor(self, r: float) -> float:
        """
        Performance degradation due to plume overlap
        
        Args:
            r: Distance between emitters (m)
        
        Returns:
            overlap_factor: Multiplicative factor (0-1)
        
        Physics:
        - Plumes diverge at angle θ
        - At distance z downstream, plume radius ~ z·tan(θ)
        - Overlap creates backpressure, reduces emission
        """
        # Characteristic plume distance
        z_plume = 0.01  # 1 cm downstream
        r_plume = z_plume * np.tan(self.divergence_angle)
        
        # Overlap fraction
        overlap = max(0, 1 - r / (2 * r_plume))
        
        # Performance reduction (empirical: 30% max degradation at full overlap)
        factor = 1.0 - 0.3 * overlap
        
        return max(factor, 0.7)  # Never worse than 70%
    
    def thermal_coupling_resistance(self, r: float, substrate_k: float = 150) -> float:
        """
        Thermal resistance between emitters through substrate
        
        Args:
            r: Distance (m)
            substrate_k: Thermal conductivity (W/(m·K)) - default Silicon
        
        Returns:
            R_thermal: Thermal resistance (K/W)
        
        Physics:
        R = 1/(4πk·r) for spreading resistance in semi-infinite medium
        """
        if r < 1e-10:
            return 1e10  # Infinite resistance (same point)
        
        R = 1 / (4 * np.pi * substrate_k * r)
        
        return R
    
    def radiative_coupling(self, r: float, T1: float, T2: float,
                          emissivity: float = 0.8) -> float:
        """
        Radiative heat transfer between emitters
        
        Args:
            r: Distance (m)
            T1, T2: Temperatures (K)
            emissivity: Surface emissivity
        
        Returns:
            Q_rad: Heat flow from 1 to 2 (W)
        
        Physics:
        Q = σ·ε·A·F12·(T1⁴ - T2⁴)
        where F12 is view factor
        """
        # Emitter area (approximate)
        A = np.pi * (self.d_tip / 2)**2
        
        # View factor (for small disks at distance r)
        # Approximate: F ~ A/πr² for r >> sqrt(A)
        F12 = A / (np.pi * r**2) if r > 10 * np.sqrt(A) else 0.01
        
        # Radiative heat transfer
        Q = PhysicalConstants.sigma_SB * emissivity * A * F12 * (T1**4 - T2**4)
        
        return Q
    
    def calculate_total_coupling_factor(self, emitter_i: int,
                                       positions: np.ndarray,
                                       firing_pattern: np.ndarray,
                                       temperatures: Optional[np.ndarray] = None) -> float:
        """
        Calculate total performance factor for emitter i
        accounting for ALL neighbors
        
        Args:
            emitter_i: Index of emitter to analyze
            positions: Nx3 array of all positions
            firing_pattern: Boolean array of which are active
            temperatures: Optional temperature array
        
        Returns:
            total_factor: Combined coupling factor (0-1)
        """
        pos_i = positions[emitter_i]
        
        # Initialize factors
        em_factor = 1.0
        plume_factor = 1.0
        
        # Loop over all other active emitters
        for j, active in enumerate(firing_pattern):
            if j == emitter_i or not active:
                continue
            
            pos_j = positions[j]
            r_vec = pos_i - pos_j
            r = np.linalg.norm(r_vec)
            
            if r < 1e-10:
                continue
            
            # EM coupling
            em_factor *= self.electric_field_coupling(r)
            
            # Plume overlap
            plume_factor *= self.plume_overlap_factor(r)
        
        # Combined effect (multiplicative)
        total_factor = em_factor * plume_factor
        
        # Thermal effects (simplified - would need iteration for full coupling)
        if temperatures is not None:
            T_i = temperatures[emitter_i]
            # Conductivity enhancement with temperature
            # K(T) ~ exp(-Ea/(RT))
            # Approximate: 5% increase per 100K
            temp_factor = 1.0 + 0.0005 * (T_i - 300)
            total_factor *= temp_factor
        
        return total_factor


# ============================================================================
# ARRAY PERFORMANCE CALCULATOR
# ============================================================================

class ArrayPerformanceCalculator:
    """
    Calculate total array performance including all coupling effects
    """
    
    def __init__(self, coupling_physics: EmitterCouplingPhysics):
        self.physics = coupling_physics
    
    def calculate_performance(self, positions: np.ndarray,
                             firing_pattern: Optional[np.ndarray] = None) -> Dict:
        """
        Calculate complete array performance
        
        Args:
            positions: Nx3 emitter positions
            firing_pattern: Which emitters fire (default: all)
        
        Returns:
            performance: Dict with all metrics
        """
        n = len(positions)
        
        if firing_pattern is None:
            firing_pattern = np.ones(n, dtype=bool)
        
        # Calculate temperatures (simplified - no full thermal network)
        # Assume each emitter at operating temp + coupling
        T = 300 + 100 * firing_pattern  # Base: 300K + 100K rise when on
        
        # Calculate individual emitter performance with coupling
        thrust_individual = np.zeros(n)
        current_individual = np.zeros(n)
        
        for i in range(n):
            if firing_pattern[i]:
                # Coupling factor from all neighbors
                factor = self.physics.calculate_total_coupling_factor(
                    i, positions, firing_pattern, T
                )
                
                thrust_individual[i] = self.physics.thrust_single * factor
                current_individual[i] = self.physics.I_single * factor
        
        # Totals
        thrust_total = np.sum(thrust_individual)
        current_total = np.sum(current_individual)
        power_total = current_total * self.physics.V
        
        # Array metrics
        array_diameter = 2 * np.max(np.linalg.norm(positions[:, :2], axis=1))
        array_area = np.pi * (array_diameter / 2)**2
        
        thrust_density = thrust_total / array_area if array_area > 0 else 0
        
        # Average coupling factor
        avg_factor = np.mean([
            self.physics.calculate_total_coupling_factor(i, positions, firing_pattern, T)
            for i in range(n) if firing_pattern[i]
        ]) if np.any(firing_pattern) else 0
        
        return {
            'thrust_total': thrust_total,
            'current_total': current_total,
            'power_total': power_total,
            'thrust_individual': thrust_individual,
            'n_active': np.sum(firing_pattern),
            'array_diameter': array_diameter,
            'array_area': array_area,
            'thrust_density': thrust_density,
            'avg_coupling_factor': avg_factor,
            'efficiency_vs_ideal': avg_factor  # How close to ideal (no coupling)
        }


# ============================================================================
# GEOMETRY OPTIMIZATION
# ============================================================================

class GeometryOptimizer:
    """
    Optimize emitter arrangement - NOT LIMITED TO RECTANGLES!
    
    Uses global optimization to find best positions
    """
    
    def __init__(self, n_emitters: int, calculator: ArrayPerformanceCalculator,
                 constraints: Dict):
        """
        Args:
            n_emitters: Number of emitters to place
            calculator: Performance calculator
            constraints: Dict with max_diameter, min_spacing, etc.
        """
        self.n_emitters = n_emitters
        self.calculator = calculator
        self.constraints = constraints
        
        self.max_diameter = constraints.get('max_diameter', 50e-3)
        self.min_spacing = constraints.get('min_spacing', 0.5e-3)
        
        print(f"\nGeometry Optimizer initialized:")
        print(f"  Emitters: {n_emitters}")
        print(f"  Max diameter: {self.max_diameter*1e3:.1f} mm")
        print(f"  Min spacing: {self.min_spacing*1e3:.2f} mm")
    
    def generate_initial_guess(self, config_type: str = 'hexagonal') -> np.ndarray:
        """
        Generate initial guess for optimization
        
        Args:
            config_type: 'hexagonal', 'square', 'circular', 'random'
        """
        if config_type == 'hexagonal':
            return self._hexagonal_grid()
        elif config_type == 'square':
            return self._square_grid()
        elif config_type == 'circular':
            return self._circular_grid()
        elif config_type == 'random':
            return self._random_grid()
        else:
            return self._hexagonal_grid()
    
    def _hexagonal_grid(self, spacing: float = 2e-3) -> np.ndarray:
        """Hexagonal close-packed grid"""
        positions = [[0, 0, 0]]  # Center
        
        # Add rings
        ring = 1
        while len(positions) < self.n_emitters:
            for i in range(6 * ring):
                if len(positions) >= self.n_emitters:
                    break
                
                # Hexagonal coordinates
                angle_step = 2 * np.pi / (6 * ring)
                angle = i * angle_step
                
                x = ring * spacing * np.cos(angle)
                y = ring * spacing * np.sin(angle)
                
                positions.append([x, y, 0])
            
            ring += 1
        
        return np.array(positions[:self.n_emitters])
    
    def _square_grid(self, spacing: float = 2e-3) -> np.ndarray:
        """Square grid"""
        n_side = int(np.ceil(np.sqrt(self.n_emitters)))
        positions = []
        
        for i in range(n_side):
            for j in range(n_side):
                if len(positions) >= self.n_emitters:
                    break
                
                x = (i - n_side/2 + 0.5) * spacing
                y = (j - n_side/2 + 0.5) * spacing
                
                positions.append([x, y, 0])
        
        return np.array(positions[:self.n_emitters])
    
    def _circular_grid(self, radius: float = 5e-3) -> np.ndarray:
        """Circular rings"""
        positions = [[0, 0, 0]]
        
        n_rings = int(np.sqrt(self.n_emitters / 3))
        
        for ring in range(1, n_rings + 1):
            r = radius * ring / n_rings
            n_in_ring = min(6 * ring, self.n_emitters - len(positions))
            
            for i in range(n_in_ring):
                angle = 2 * np.pi * i / n_in_ring
                x = r * np.cos(angle)
                y = r * np.sin(angle)
                positions.append([x, y, 0])
                
                if len(positions) >= self.n_emitters:
                    break
        
        return np.array(positions[:self.n_emitters])
    
    def _random_grid(self) -> np.ndarray:
        """Random positions (for optimization starting point)"""
        r_max = self.max_diameter / 2 * 0.8
        
        positions = []
        for i in range(self.n_emitters):
            r = r_max * np.random.random()
            theta = 2 * np.pi * np.random.random()
            
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            
            positions.append([x, y, 0])
        
        return np.array(positions)
    
    def objective_function(self, x: np.ndarray) -> float:
        """
        Objective function for optimization
        
        Args:
            x: Flattened array of positions [x1,y1,x2,y2,...,xN,yN]
        
        Returns:
            objective: Value to minimize (negative performance)
        """
        # Reshape to positions
        positions = x.reshape(-1, 2)
        positions_3d = np.column_stack([positions, np.zeros(len(positions))])
        
        # Check constraints
        # 1. Minimum spacing
        dist_matrix = distance_matrix(positions, positions)
        np.fill_diagonal(dist_matrix, np.inf)
        min_dist = np.min(dist_matrix)
        
        if min_dist < self.min_spacing:
            return 1e10  # Penalty
        
        # 2. Maximum diameter
        max_r = np.max(np.linalg.norm(positions, axis=1))
        diameter = 2 * max_r
        
        if diameter > self.max_diameter:
            return 1e10  # Penalty
        
        # Calculate performance
        perf = self.calculator.calculate_performance(positions_3d)
        
        # Maximize thrust (minimize negative thrust)
        objective = -perf['thrust_total']
        
        # Bonus for compact arrays (thrust density)
        # objective += 1e-9 / perf['thrust_density'] if perf['thrust_density'] > 0 else 0
        
        return objective
    
    def optimize(self, method: str = 'differential_evolution',
                initial_config: str = 'hexagonal') -> Tuple[np.ndarray, Dict]:
        """
        Run optimization
        
        Args:
            method: 'differential_evolution', 'nelder-mead', etc.
            initial_config: Initial guess type
        
        Returns:
            (optimal_positions, performance)
        """
        print(f"\n{'='*70}")
        print(f"RUNNING GEOMETRY OPTIMIZATION")
        print(f"{'='*70}")
        print(f"Method: {method}")
        print(f"Initial config: {initial_config}")
        
        # Initial guess
        x0 = self.generate_initial_guess(initial_config)
        x0_flat = x0[:, :2].flatten()
        
        # Evaluate initial
        perf0 = self.calculator.calculate_performance(x0)
        print(f"\nInitial performance:")
        print(f"  Thrust: {perf0['thrust_total']*1e6:.2f} μN")
        print(f"  Diameter: {perf0['array_diameter']*1e3:.2f} mm")
        print(f"  Coupling factor: {perf0['avg_coupling_factor']:.3f}")
        
        # Bounds (positions within max_diameter/2)
        r_max = self.max_diameter / 2
        bounds = [(-r_max, r_max)] * (2 * self.n_emitters)
        
        # Optimize
        print(f"\nOptimizing...")
        
        if method == 'differential_evolution':
            result = differential_evolution(
                self.objective_function,
                bounds=bounds,
                maxiter=100,
                popsize=15,
                tol=1e-7,
                seed=42,
                workers=1,
                updating='deferred',
                polish=True
            )
        else:
            result = minimize(
                self.objective_function,
                x0=x0_flat,
                method='Nelder-Mead',
                options={'maxiter': 1000, 'xatol': 1e-6}
            )
        
        # Extract result
        x_opt = result.x.reshape(-1, 2)
        positions_opt = np.column_stack([x_opt, np.zeros(len(x_opt))])
        
        # Final performance
        perf_opt = self.calculator.calculate_performance(positions_opt)
        
        print(f"\nOptimization complete!")
        print(f"  Status: {result.message if hasattr(result, 'message') else 'success'}")
        print(f"  Iterations: {result.nit if hasattr(result, 'nit') else 'N/A'}")
        print(f"\nOptimal performance:")
        print(f"  Thrust: {perf_opt['thrust_total']*1e6:.2f} μN "
              f"({(perf_opt['thrust_total']/perf0['thrust_total']-1)*100:+.1f}% vs initial)")
        print(f"  Diameter: {perf_opt['array_diameter']*1e3:.2f} mm")
        print(f"  Coupling factor: {perf_opt['avg_coupling_factor']:.3f}")
        print(f"  Thrust density: {perf_opt['thrust_density']*1e9:.2f} μN/mm²")
        
        return positions_opt, perf_opt


# ============================================================================
# TO BE CONTINUED IN PART 2
# ============================================================================
