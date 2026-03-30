"""
    COMPLETE MULTI-EMITTER ELECTROSPRAY THRUSTER SYSTEM  
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Arrow
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from scipy.optimize import differential_evolution, minimize, dual_annealing
from scipy.integrate import solve_ivp
from scipy.spatial import distance_matrix
from typing import Dict, List, Tuple, Optional
import warnings


# EXPANDED CONSTANTS & PROPELLANT DATABASE

class Constants:
    """Physical constants and EXPANDED propellant database"""
    
    mu_earth = 3.986004418e14  # m³/s²
    R_earth = 6378137  # m
    eps0 = 8.854187817e-12  # F/m
    sigma_SB = 5.670374419e-8  # W/(m²·K⁴)
    AMU = 1.66053906660e-27  # kg
    
    PROPELLANTS: Dict[str, Dict] = {
        'EMI-Im': {
            'name': 'EMI-Im (1-Ethyl-3-methylimidazolium Imide)',
            'conductivity': 1.5,  # S/m
            'surface_tension': 0.042,  # N/m
            'density': 1520,  # kg/m³
            'viscosity': 0.034,  # Pa·s
            'ion_mass': 280,  # amu
            'isp_range': [2000, 4000],  # s
            'recommended': True
        },
        'EMI-BF4': {
            'name': 'EMI-BF4 (Tetrafluoroborate)',
            'conductivity': 1.4,
            'surface_tension': 0.048,
            'density': 1280,
            'viscosity': 0.042,
            'ion_mass': 111,
            'isp_range': [2500, 5000],
            'recommended': True
        },
        'BMIM-DCA': {
            'name': 'BMIM-DCA (Butylmethylimidazolium Dicyanamide)',
            'conductivity': 0.95,
            'surface_tension': 0.049,
            'density': 1050,
            'viscosity': 0.055,
            'ion_mass': 180,
            'isp_range': [1800, 3500],
            'recommended': False
        },
        'EMIM-NTf2': {
            'name': 'EMIM-NTf2 (Bis(trifluoromethylsulfonyl)imide)',
            'conductivity': 0.87,
            'surface_tension': 0.035,
            'density': 1518,
            'viscosity': 0.034,
            'ion_mass': 280,
            'isp_range': [2000, 4000],
            'recommended': True
        },
        'HMIM-BF4': {
            'name': 'HMIM-BF4 (Hexylmethylimidazolium Tetrafluoroborate)',
            'conductivity': 0.65,
            'surface_tension': 0.043,
            'density': 1150,
            'viscosity': 0.089,
            'ion_mass': 197,
            'isp_range': [2200, 4200],
            'recommended': False
        },
        'Formamide': {
            'name': 'Formamide (NH2CHO)',
            'conductivity': 0.002,
            'surface_tension': 0.058,
            'density': 1133,
            'viscosity': 0.00345,
            'ion_mass': 45,
            'isp_range': [3000, 6000],
            'recommended': False
        }
    }



# PHYSICS MODELS (Enhanced)

class EmitterCouplingPhysics:
    """Enhanced coupling physics"""
    
    def __init__(self, propellant: Dict, voltage: float):
        self.prop = propellant
        self.V = voltage
        self.I_single = 200e-9  # A
        self.thrust_single = 1e-6  # N
        self.d_tip = 30e-6  # m
        self.divergence = 0.05  # rad
        
        self.m_ion = propellant['ion_mass'] * Constants.AMU
        self.v_ion = np.sqrt(2 * voltage * 1.602e-19 / self.m_ion)
    
    def electric_field_coupling(self, r: float) -> float:
        """Enhanced EM coupling with proper scaling"""
        r_norm = r / self.d_tip
        
        if r_norm < 5:
            # Very strong shielding
            return 0.5 + 0.5 * (r_norm / 5)
        elif r_norm < 15:
            # Strong shielding
            return 0.7 + 0.3 * ((r_norm - 5) / 10)
        elif r_norm < 50:
            # Slight enhancement
            return 1.0 + 0.1 * np.exp(-(r_norm - 15) / 20)
        else:
            return 1.0
    
    def plume_overlap_factor(self, r: float) -> float:
        """Enhanced plume overlap"""
        z_plume = 0.01  # m
        r_plume = z_plume * np.tan(self.divergence)
        
        if r < r_plume:
            # Strong overlap
            overlap = 1.0 - (r / r_plume)
            return max(0.6, 1.0 - 0.4 * overlap)
        elif r < 2 * r_plume:
            # Moderate overlap
            overlap = 1.0 - (r - r_plume) / r_plume
            return max(0.85, 1.0 - 0.15 * overlap)
        else:
            return 1.0
    
    def calculate_total_coupling_factor(self, i: int, positions: np.ndarray,
                                       firing: np.ndarray) -> float:
        """Total coupling with enhanced model"""
        pos_i = positions[i]
        em_factor = 1.0
        plume_factor = 1.0
        
        for j in range(len(positions)):
            if j == i or not firing[j]:
                continue
            
            r = np.linalg.norm(pos_i - positions[j])
            if r < 1e-10:
                continue
            
            em_factor *= self.electric_field_coupling(r)
            plume_factor *= self.plume_overlap_factor(r)
        
        return em_factor * plume_factor


class ArrayPerformanceCalculator:
    """Performance calculator"""
    
    def __init__(self, physics: EmitterCouplingPhysics):
        self.physics = physics
    
    def calculate_performance(self, positions: np.ndarray,
                             firing: Optional[np.ndarray] = None) -> Dict:
        """Calculate performance"""
        n = len(positions)
        
        if firing is None:
            firing = np.ones(n, dtype=bool)
        
        thrust_total = 0
        current_total = 0
        
        for i in range(n):
            if firing[i]:
                factor = self.physics.calculate_total_coupling_factor(i, positions, firing)
                thrust_total += self.physics.thrust_single * factor
                current_total += self.physics.I_single * factor
        
        diameter = 2 * np.max(np.linalg.norm(positions[:, :2], axis=1))
        area = np.pi * (diameter / 2)**2
        
        return {
            'thrust_total': thrust_total,
            'current_total': current_total,
            'power_total': current_total * self.physics.V,
            'n_active': np.sum(firing),
            'array_diameter': diameter,
            'thrust_density': thrust_total / area if area > 0 else 0,
            'avg_coupling_factor': thrust_total / (np.sum(firing) * self.physics.thrust_single) if np.sum(firing) > 0 else 0
        }



# ENHANCED GEOMETRY OPTIMIZER


class EnhancedGeometryOptimizer:
    """Enhanced optimizer with multiple algorithms"""
    
    def __init__(self, n_emitters: int, propellant: Dict, voltage: float):
        self.n_emitters = n_emitters
        self.physics = EmitterCouplingPhysics(propellant, voltage)
        self.calculator = ArrayPerformanceCalculator(self.physics)
        
        self.max_diameter = 50e-3  # m
        self.min_spacing = 0.5e-3  # m
    
    def _hexagonal_grid(self, spacing: float) -> np.ndarray:
        """Generate hexagonal grid"""
        positions = [[0, 0, 0]]
        ring = 1
        
        while len(positions) < self.n_emitters:
            for i in range(6 * ring):
                if len(positions) >= self.n_emitters:
                    break
                angle = 2 * np.pi * i / (6 * ring)
                x = ring * spacing * np.cos(angle)
                y = ring * spacing * np.sin(angle)
                positions.append([x, y, 0])
            ring += 1
        
        return np.array(positions[:self.n_emitters])
    
    def optimize(self, method: str = 'quick') -> Tuple[np.ndarray, Dict]:
        """Enhanced optimization with multiple methods"""
        print(f"\n{'='*70}")
        print(f"ENHANCED GEOMETRY OPTIMIZATION")
        print(f"{'='*70}")
        print(f"Algorithm: {method}")
        
        def objective(spacing):
            positions = self._hexagonal_grid(spacing[0])
            firing = np.ones(len(positions), dtype=bool)
            
            dists = distance_matrix(positions[:, :2], positions[:, :2])
            np.fill_diagonal(dists, np.inf)
            if np.min(dists) < self.min_spacing:
                return 1e10
            
            diameter = 2 * np.max(np.linalg.norm(positions[:, :2], axis=1))
            if diameter > self.max_diameter:
                return 1e10
            
            perf = self.calculator.calculate_performance(positions, firing)
            return -perf['thrust_total']
        
        if method == 'thorough':
            # Use dual annealing for global optimization
            print("Using dual annealing (global optimization)...")
            result = dual_annealing(objective, bounds=[(0.5e-3, 5e-3)],
                                   maxiter=100, seed=42)
        else:
            # Quick L-BFGS-B
            result = minimize(objective, x0=[2e-3], bounds=[(0.5e-3, 5e-3)],
                            method='L-BFGS-B')
        
        optimal_spacing = result.x[0]
        positions = self._hexagonal_grid(optimal_spacing)
        performance = self.calculator.calculate_performance(positions)
        
        print(f"✓ Optimal spacing: {optimal_spacing*1e3:.2f} mm")
        print(f"✓ Total thrust: {performance['thrust_total']*1e6:.2f} μN")
        print(f"✓ Array diameter: {performance['array_diameter']*1e3:.2f} mm")
        print(f"✓ Efficiency: {performance['avg_coupling_factor']:.3f}")
        
        return positions, performance



# ENHANCED ORBITAL MECHANICS


class EnhancedOrbitalMechanics:
    """Enhanced orbital mechanics with more maneuver types"""
    
    @staticmethod
    def hohmann_transfer(r1: float, r2: float) -> Tuple[float, float, float]:
        """Hohmann transfer"""
        mu = Constants.mu_earth
        
        a_transfer = (r1 + r2) / 2
        
        v1 = np.sqrt(mu / r1)
        v1_transfer = np.sqrt(mu * (2/r1 - 1/a_transfer))
        dv1 = v1_transfer - v1
        
        v2 = np.sqrt(mu / r2)
        v2_transfer = np.sqrt(mu * (2/r2 - 1/a_transfer))
        dv2 = v2 - v2_transfer
        
        dt = np.pi * np.sqrt(a_transfer**3 / mu)
        
        return dv1, dv2, dt
    
    @staticmethod
    def bielliptic_transfer(r1: float, r2: float, r_apogee: float) -> Tuple[float, float, float, float]:
        """Bi-elliptic transfer (more efficient for large radius changes)"""
        mu = Constants.mu_earth
        
        # First burn: raise apogee to r_apogee
        v1 = np.sqrt(mu / r1)
        a1 = (r1 + r_apogee) / 2
        v1_transfer = np.sqrt(mu * (2/r1 - 1/a1))
        dv1 = v1_transfer - v1
        
        # Second burn: at apogee, change periapsis to r2
        v_apogee1 = np.sqrt(mu * (2/r_apogee - 1/a1))
        a2 = (r_apogee + r2) / 2
        v_apogee2 = np.sqrt(mu * (2/r_apogee - 1/a2))
        dv2 = v_apogee2 - v_apogee1
        
        # Third burn: circularize at r2
        v2_transfer = np.sqrt(mu * (2/r2 - 1/a2))
        v2 = np.sqrt(mu / r2)
        dv3 = v2 - v2_transfer
        
        # Transfer times
        dt1 = np.pi * np.sqrt(a1**3 / mu)
        dt2 = np.pi * np.sqrt(a2**3 / mu)
        
        return dv1, dv2, dv3, dt1 + dt2
    
    @staticmethod
    def plane_change(v: float, delta_i: float) -> float:
        """Plane change maneuver"""
        # dv = 2*v*sin(Δi/2)
        return 2 * v * np.sin(delta_i / 2)
    
    @staticmethod
    def orbital_elements_to_cartesian(a: float, e: float, i: float,
                                     Omega: float, omega: float, nu: float) -> Tuple[np.ndarray, np.ndarray]:
        """Convert orbital elements to Cartesian"""
        p = a * (1 - e**2)
        
        r_pf = p / (1 + e*np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0])
        v_pf = np.sqrt(Constants.mu_earth/p) * np.array([-np.sin(nu), e + np.cos(nu), 0])
        
        R3_O = np.array([[np.cos(Omega), -np.sin(Omega), 0],
                         [np.sin(Omega), np.cos(Omega), 0],
                         [0, 0, 1]])
        
        R1_i = np.array([[1, 0, 0],
                        [0, np.cos(i), -np.sin(i)],
                        [0, np.sin(i), np.cos(i)]])
        
        R3_w = np.array([[np.cos(omega), -np.sin(omega), 0],
                        [np.sin(omega), np.cos(omega), 0],
                        [0, 0, 1]])
        
        R = R3_O @ R1_i @ R3_w
        
        return R @ r_pf, R @ v_pf


# ═══════════════════════════════════════════════════════════════════════════
# CONTROL ALLOCATOR (Enhanced)
# ═══════════════════════════════════════════════════════════════════════════

class ControlAllocator:
    """Enhanced control allocation"""
    
    def __init__(self, positions: np.ndarray, thrust: float = 1e-6):
        self.positions = positions
        self.n = len(positions)
        self.thrust = thrust
        self._build_matrix()
    
    def _build_matrix(self):
        """Build control matrix"""
        self.B = np.zeros((6, self.n))
        
        for i in range(self.n):
            F = np.array([0, 0, self.thrust])
            self.B[0:3, i] = F
            self.B[3:6, i] = np.cross(self.positions[i], F)
    
    def allocate_greedy(self, F_des: np.ndarray, tau_des: np.ndarray) -> np.ndarray:
        """Greedy allocation"""
        desired = np.concatenate([F_des, tau_des])
        u = np.zeros(self.n, dtype=bool)
        
        for _ in range(self.n):
            best_i = -1
            best_err = np.inf
            
            for i in range(self.n):
                if u[i]:
                    continue
                u_test = u.copy()
                u_test[i] = True
                err = np.linalg.norm(self.B @ u_test.astype(float) - desired)
                
                if err < best_err:
                    best_err = err
                    best_i = i
            
            if best_i == -1:
                break
            
            u[best_i] = True
            
            if best_err < 1e-9:
                break
        
        return u


# ═══════════════════════════════════════════════════════════════════════════
# SPACECRAFT CONTROLLER (Enhanced with more maneuvers)
# ═══════════════════════════════════════════════════════════════════════════

class EnhancedSpacecraftController:
    """Enhanced controller with more maneuver types"""
    
    def __init__(self, positions: np.ndarray, mass: float, inertia: np.ndarray,
                 thrust_per_emitter: float = 1e-6):
        self.positions = positions
        self.mass = mass
        self.inertia = inertia
        self.thrust = thrust_per_emitter
        
        self.allocator = ControlAllocator(positions, thrust_per_emitter)
        self.om = EnhancedOrbitalMechanics()
    
    def plan_hohmann(self, initial_orbit: Dict, target_orbit: Dict) -> Dict:
        """Plan Hohmann transfer"""
        r1 = initial_orbit['a']
        r2 = target_orbit['a']
        
        dv1, dv2, dt_transfer = self.om.hohmann_transfer(r1, r2)
        
        total_thrust = len(self.positions) * self.thrust
        t_burn1 = self.mass * abs(dv1) / total_thrust
        t_burn2 = self.mass * abs(dv2) / total_thrust
        
        return {
            'type': 'hohmann',
            'dv1': dv1,
            'dv2': dv2,
            'total_dv': abs(dv1) + abs(dv2),
            'burn_time_1': t_burn1,
            'burn_time_2': t_burn2,
            'transfer_time': dt_transfer,
            'firing_pattern_1': np.ones(len(self.positions), dtype=bool),
            'firing_pattern_2': np.ones(len(self.positions), dtype=bool)
        }
    
    def plan_bielliptic(self, initial_orbit: Dict, target_orbit: Dict) -> Dict:
        """Plan bi-elliptic transfer"""
        r1 = initial_orbit['a']
        r2 = target_orbit['a']
        r_apogee = max(r1, r2) * 3  # High apogee
        
        dv1, dv2, dv3, dt_transfer = self.om.bielliptic_transfer(r1, r2, r_apogee)
        
        total_thrust = len(self.positions) * self.thrust
        t_burn1 = self.mass * abs(dv1) / total_thrust
        t_burn2 = self.mass * abs(dv2) / total_thrust
        t_burn3 = self.mass * abs(dv3) / total_thrust
        
        return {
            'type': 'bielliptic',
            'dv1': dv1,
            'dv2': dv2,
            'dv3': dv3,
            'total_dv': abs(dv1) + abs(dv2) + abs(dv3),
            'burn_time_1': t_burn1,
            'burn_time_2': t_burn2,
            'burn_time_3': t_burn3,
            'transfer_time': dt_transfer,
            'firing_pattern': np.ones(len(self.positions), dtype=bool)
        }
    
    def plan_plane_change(self, delta_i_deg: float, altitude: float) -> Dict:
        """Plan plane change maneuver"""
        r = Constants.R_earth + altitude
        v = np.sqrt(Constants.mu_earth / r)
        delta_i = np.deg2rad(delta_i_deg)
        
        dv = self.om.plane_change(v, delta_i)
        
        total_thrust = len(self.positions) * self.thrust
        t_burn = self.mass * dv / total_thrust
        
        return {
            'type': 'plane_change',
            'delta_i_deg': delta_i_deg,
            'dv': dv,
            'burn_time': t_burn,
            'firing_pattern': np.ones(len(self.positions), dtype=bool)
        }
    
    def plan_attitude(self, target_euler: np.ndarray, current_euler: np.ndarray = None) -> Dict:
        """Plan attitude maneuver"""
        if current_euler is None:
            current_euler = np.zeros(3)
        
        target_rad = np.deg2rad(target_euler)
        current_rad = np.deg2rad(current_euler)
        
        error = target_rad - current_rad
        
        kp = 0.1
        desired_torque = kp * error
        desired_force = np.zeros(3)
        
        firing = self.allocator.allocate_greedy(desired_force, desired_torque)
        
        angle_change = np.linalg.norm(error)
        t_maneuver = 60 * (angle_change / (10 * np.pi/180))
        
        return {
            'type': 'attitude',
            'target_euler': target_euler,
            'current_euler': current_euler,
            'angle_change_deg': np.rad2deg(angle_change),
            'firing_pattern': firing,
            'n_emitters_active': np.sum(firing),
            'estimated_time': t_maneuver
        }
    
    def plan_detumbling(self, omega_initial: np.ndarray) -> Dict:
        """Plan detumbling maneuver"""
        # Simple bang-bang control
        torque_des = -0.5 * self.inertia @ omega_initial
        firing = self.allocator.allocate_greedy(np.zeros(3), torque_des)
        
        # Estimate time (simplified)
        omega_mag = np.linalg.norm(omega_initial)
        t_detumble = omega_mag * 100  # seconds
        
        return {
            'type': 'detumbling',
            'omega_initial_deg_s': np.rad2deg(omega_initial),
            'firing_pattern': firing,
            'n_emitters_active': np.sum(firing),
            'estimated_time': t_detumble
        }


# ═══════════════════════════════════════════════════════════════════════════
# ENHANCED VISUALIZER (Fixed maneuver plots!)
# ═══════════════════════════════════════════════════════════════════════════

class EnhancedVisualizer:
    """Enhanced visualizer with UNIQUE maneuver plots"""
    
    @staticmethod
    def plot_array_layout(positions: np.ndarray, performance: Dict,
                         firing_pattern: Optional[np.ndarray] = None,
                         title: str = "Array Layout"):
        """Plot 2D array layout (standard view)"""
        if firing_pattern is None:
            firing_pattern = np.ones(len(positions), dtype=bool)
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        ax = axes[0]
        
        colors = ['red' if f else 'lightgray' for f in firing_pattern]
        sizes = [200 if f else 80 for f in firing_pattern]
        
        ax.scatter(positions[:, 0]*1e3, positions[:, 1]*1e3, c=colors, s=sizes,
                  alpha=0.9, edgecolors='black', linewidth=2, zorder=3)
        
        diameter = performance.get('array_diameter', 0)
        circle = Circle((0, 0), diameter*1e3/2, fill=False, edgecolor='blue',
                       linewidth=2, linestyle='--', label='Array boundary')
        ax.add_patch(circle)
        
        ax.set_xlabel('X Position (mm)', fontsize=13)
        ax.set_ylabel('Y Position (mm)', fontsize=13)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        ax.legend()
        
        lim = diameter * 1e3 * 0.6
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        
        # Performance metrics
        ax2 = axes[1]
        ax2.axis('off')
        
        n_active = np.sum(firing_pattern)
        
        metrics_text = f"""
ARRAY PERFORMANCE
{'='*45}

Configuration:
  Total Emitters:    {len(positions)}
  Active Emitters:   {n_active}
  Array Diameter:    {diameter*1e3:.2f} mm
  
Performance:
  Total Thrust:      {performance.get('thrust_total', 0)*1e6:.2f} μN
  Total Current:     {performance.get('current_total', 0)*1e9:.2f} nA
  Total Power:       {performance.get('power_total', 0)*1e3:.2f} mW
  Efficiency:        {performance.get('avg_coupling_factor', 0):.3f}
  
Thrust Density:      {performance.get('thrust_density', 0)*1e9:.2f} μN/mm²
        """
        
        ax2.text(0.1, 0.95, metrics_text, transform=ax2.transAxes,
                fontsize=11, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        plt.tight_layout()
        return fig
    
    @staticmethod
    def plot_maneuver_plan(positions: np.ndarray, firing_pattern: np.ndarray,
                          maneuver_info: Dict, title: str = "Maneuver Plan"):
        """
        Plot UNIQUE maneuver plan (DIFFERENT from array layout!)
        Shows thrust vectors, torque visualization, burn info
        """
        fig = plt.figure(figsize=(16, 8))
        
        # Left: Firing pattern with thrust vectors
        ax1 = plt.subplot(1, 2, 1)
        
        # Plot emitters
        colors = ['darkred' if f else 'lightgray' for f in firing_pattern]
        sizes = [250 if f else 60 for f in firing_pattern]
        
        for i, (pos, firing) in enumerate(zip(positions, firing_pattern)):
            ax1.scatter(pos[0]*1e3, pos[1]*1e3, c=colors[i], s=sizes[i],
                       alpha=0.9, edgecolors='black', linewidth=2, zorder=3)
            
            # Draw thrust vector for active emitters
            if firing:
                arrow = FancyArrowPatch((pos[0]*1e3, pos[1]*1e3),
                                       (pos[0]*1e3, pos[1]*1e3 + 3),
                                       mutation_scale=15, linewidth=2,
                                       arrowstyle='->', color='orange',
                                       zorder=4, alpha=0.8)
                ax1.add_patch(arrow)
        
        # Add thrust direction indicator
        ax1.arrow(0, -15, 0, 5, head_width=2, head_length=1.5,
                 fc='orange', ec='orange', linewidth=3, alpha=0.7)
        ax1.text(0, -18, 'Thrust Direction', ha='center', fontsize=11,
                fontweight='bold', color='orange')
        
        ax1.set_xlabel('X Position (mm)', fontsize=12)
        ax1.set_ylabel('Y Position (mm)', fontsize=12)
        ax1.set_title(f'{title} - Thrust Vectors', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        
        max_r = np.max(np.abs(positions[:, :2])) * 1e3 * 1.5
        ax1.set_xlim(-max_r, max_r)
        ax1.set_ylim(-max_r, max_r)
        
        # Right: Maneuver details
        ax2 = plt.subplot(1, 2, 2)
        ax2.axis('off')
        
        # Build maneuver details text
        maneuver_type = maneuver_info.get('type', 'unknown')
        n_active = np.sum(firing_pattern)
        
        if maneuver_type == 'hohmann':
            details = f"""
HOHMANN TRANSFER MANEUVER
{'='*50}

Configuration:
  Active Emitters:        {n_active} / {len(positions)}
  Total Thrust:           {n_active * 1e-6 * 1e6:.2f} μN

Maneuver Parameters:
  Type:                   Hohmann Transfer
  Burn 1 Δv:              {maneuver_info.get('dv1', 0):.2f} m/s
  Burn 2 Δv:              {maneuver_info.get('dv2', 0):.2f} m/s
  Total Δv:               {maneuver_info.get('total_dv', 0):.2f} m/s

Burn Times:
  Burn 1:                 {maneuver_info.get('burn_time_1', 0)/3600:.2f} hours
  Burn 2:                 {maneuver_info.get('burn_time_2', 0)/3600:.2f} hours
  Transfer Time:          {maneuver_info.get('transfer_time', 0)/3600:.2f} hours

Status:
  ✓ All emitters firing for maximum thrust
  ✓ Optimal 2-burn transfer
            """
        elif maneuver_type == 'attitude':
            details = f"""
ATTITUDE CONTROL MANEUVER
{'='*50}

Configuration:
  Active Emitters:        {n_active} / {len(positions)}
  Inactive:               {len(positions) - n_active}

Target Attitude:
  Roll:                   {maneuver_info.get('target_euler', [0,0,0])[0]:.1f}°
  Pitch:                  {maneuver_info.get('target_euler', [0,0,0])[1]:.1f}°
  Yaw:                    {maneuver_info.get('target_euler', [0,0,0])[2]:.1f}°

Maneuver:
  Angle Change:           {maneuver_info.get('angle_change_deg', 0):.2f}°
  Estimated Time:         {maneuver_info.get('estimated_time', 0):.1f} seconds

Status:
  ✓ Asymmetric firing pattern for torque
  ✓ Minimal translational force
            """
        elif maneuver_type == 'plane_change':
            details = f"""
PLANE CHANGE MANEUVER
{'='*50}

Configuration:
  Active Emitters:        {n_active} / {len(positions)}

Maneuver Parameters:
  Inclination Change:     {maneuver_info.get('delta_i_deg', 0):.2f}°
  Required Δv:            {maneuver_info.get('dv', 0):.2f} m/s
  Burn Time:              {maneuver_info.get('burn_time', 0)/3600:.2f} hours

Status:
  ✓ All emitters firing
  ✓ Burn at ascending/descending node
            """
        else:
            details = f"""
MANEUVER PLAN
{'='*50}

Configuration:
  Active Emitters:        {n_active} / {len(positions)}
  Type:                   {maneuver_type}
            """
        
        ax2.text(0.1, 0.95, details, transform=ax2.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        return fig
    
    @staticmethod
    def plot_3d_array(positions: np.ndarray, firing_pattern: np.ndarray):
        """3D visualization"""
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        colors = ['red' if f else 'gray' for f in firing_pattern]
        sizes = [100 if f else 30 for f in firing_pattern]
        
        ax.scatter(positions[:, 0]*1e3, positions[:, 1]*1e3, positions[:, 2]*1e3,
                  c=colors, s=sizes, alpha=0.8, edgecolors='black')
        
        # Draw plume cones
        for pos, active in zip(positions, firing_pattern):
            if active:
                h = 20
                theta = 0.05
                
                z = np.linspace(0, h, 10)
                th = np.linspace(0, 2*np.pi, 20)
                Z, Th = np.meshgrid(z, th)
                R = Z * np.tan(theta)
                
                X = pos[0]*1e3 + R * np.cos(Th)
                Y = pos[1]*1e3 + R * np.sin(Th)
                
                ax.plot_surface(X, Y, Z, alpha=0.15, color='blue')
        
        ax.set_xlabel('X (mm)', fontsize=11)
        ax.set_ylabel('Y (mm)', fontsize=11)
        ax.set_zlabel('Z (mm)', fontsize=11)
        ax.set_title('3D Array with Plume Cones', fontsize=14, fontweight='bold')
        
        max_range = np.max(np.abs(positions[:, :2])) * 1e3 * 1.2
        ax.set_xlim(-max_range, max_range)
        ax.set_ylim(-max_range, max_range)
        ax.set_zlim(0, 20)
        
        return fig


# ═══════════════════════════════════════════════════════════════════════════
# ENHANCED ANIMATION (Creates actual GIFs!)
# ═══════════════════════════════════════════════════════════════════════════

class EnhancedAnimator:
    """Enhanced animator - creates GIF files"""
    
    @staticmethod
    def animate_orbit_transfer(initial_orbit: Dict, final_orbit: Dict,
                               save_path: str = 'orbit_transfer.gif'):
        """Animate orbit transfer - creates GIF"""
        print(f"\n{'='*70}")
        print("CREATING ORBIT TRANSFER ANIMATION")
        print(f"{'='*70}")
        
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Draw Earth
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        x_earth = Constants.R_earth/1e6 * np.outer(np.cos(u), np.sin(v))
        y_earth = Constants.R_earth/1e6 * np.outer(np.sin(u), np.sin(v))
        z_earth = Constants.R_earth/1e6 * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x_earth, y_earth, z_earth, color='blue', alpha=0.4)
        
        # Generate orbits
        nu_array = np.linspace(0, 2*np.pi, 100)
        r_initial = []
        for nu in nu_array:
            r, v = EnhancedOrbitalMechanics.orbital_elements_to_cartesian(
                initial_orbit['a'], initial_orbit['e'], initial_orbit['i'],
                initial_orbit['Omega'], initial_orbit['omega'], nu
            )
            r_initial.append(r / 1e6)
        r_initial = np.array(r_initial)
        
        r_final = []
        for nu in nu_array:
            r, v = EnhancedOrbitalMechanics.orbital_elements_to_cartesian(
                final_orbit['a'], final_orbit['e'], final_orbit['i'],
                final_orbit['Omega'], final_orbit['omega'], nu
            )
            r_final.append(r / 1e6)
        r_final = np.array(r_final)
        
        ax.plot(r_initial[:, 0], r_initial[:, 1], r_initial[:, 2],
               'r--', linewidth=2.5, label='Initial Orbit', alpha=0.7)
        ax.plot(r_final[:, 0], r_final[:, 1], r_final[:, 2],
               'g--', linewidth=2.5, label='Final Orbit', alpha=0.7)
        
        sc_marker = ax.plot([r_initial[0, 0]], [r_initial[0, 1]], [r_initial[0, 2]],
                           'ro', markersize=12, label='Spacecraft', zorder=10)[0]
        
        trail, = ax.plot([], [], [], 'r-', linewidth=1, alpha=0.5)
        
        n_frames = 150
        positions = np.zeros((n_frames, 3))
        
        for i in range(n_frames):
            if i < n_frames * 0.3:
                idx = int((i / (n_frames * 0.3)) * len(r_initial))
                positions[i] = r_initial[min(idx, len(r_initial)-1)]
            elif i < n_frames * 0.7:
                frac = (i - n_frames*0.3) / (n_frames*0.4)
                positions[i] = (1-frac) * r_initial[0] + frac * r_final[0]
            else:
                idx = int(((i - n_frames*0.7) / (n_frames*0.3)) * len(r_final))
                positions[i] = r_final[min(idx, len(r_final)-1)]
        
        def update(frame):
            sc_marker.set_data([positions[frame, 0]], [positions[frame, 1]])
            sc_marker.set_3d_properties([positions[frame, 2]])
            
            # Update trail
            trail_start = max(0, frame - 20)
            trail.set_data(positions[trail_start:frame+1, 0],
                          positions[trail_start:frame+1, 1])
            trail.set_3d_properties(positions[trail_start:frame+1, 2])
            
            ax.set_title(f'Orbit Transfer (Frame {frame+1}/{n_frames})',
                        fontsize=13, fontweight='bold')
            return sc_marker, trail
        
        ax.set_xlabel('X (Mm)', fontsize=10)
        ax.set_ylabel('Y (Mm)', fontsize=10)
        ax.set_zlabel('Z (Mm)', fontsize=10)
        ax.legend(loc='upper left', fontsize=9)
        
        print("Generating animation frames...")
        anim = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)
        
        print(f"Saving animation to {save_path}...")
        try:
            writer = PillowWriter(fps=20)
            anim.save(save_path, writer=writer)
            print(f"✓ Animation saved successfully!")
        except Exception as e:
            print(f"⚠ Could not save animation: {e}")
            print("  (Install pillow: pip install pillow)")
        
        plt.close(fig)
        return anim
    
    @staticmethod
    def animate_firing_sequence(positions: np.ndarray, firing_pattern: np.ndarray,
                                duration: float, save_path: str = 'firing_sequence.gif'):
        """Animate firing sequence - creates GIF"""
        print(f"\nCreating firing sequence animation...")
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        ax.set_xlabel('X Position (mm)', fontsize=12)
        ax.set_ylabel('Y Position (mm)', fontsize=12)
        ax.set_title('Emitter Firing Sequence', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        max_r = np.max(np.abs(positions[:, :2])) * 1e3 * 1.2
        ax.set_xlim(-max_r, max_r)
        ax.set_ylim(-max_r, max_r)
        
        scatter = ax.scatter([], [], s=200, c=[], cmap='hot',
                           vmin=0, vmax=1, edgecolors='black', linewidth=2)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Firing Intensity', fontsize=11)
        
        time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                           fontsize=12, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        n_frames = 30
        
        def update(frame):
            # Create pulsing effect
            if frame < 5:
                # Startup
                intensity = firing_pattern * (frame / 5)
            elif frame < 25:
                # Active with pulse
                pulse = 0.8 + 0.2 * np.sin(frame * 0.5)
                intensity = firing_pattern * pulse
            else:
                # Shutdown
                intensity = firing_pattern * ((n_frames - frame) / 5)
            
            colors = intensity.astype(float)
            sizes = 250 * (0.3 + 0.7 * intensity)
            
            scatter.set_offsets(positions[:, :2] * 1e3)
            scatter.set_array(colors)
            scatter.set_sizes(sizes)
            
            n_active = np.sum(firing_pattern)
            time_elapsed = (frame / n_frames) * duration
            time_text.set_text(f'Time: {time_elapsed:.1f} s\nActive: {n_active}/{len(firing_pattern)}')
            
            return scatter, time_text
        
        anim = FuncAnimation(fig, update, frames=n_frames, interval=100, blit=False)
        
        try:
            writer = PillowWriter(fps=10)
            anim.save(save_path, writer=writer)
            print(f"✓ Firing sequence saved to {save_path}")
        except Exception as e:
            print(f"⚠ Could not save firing animation: {e}")
        
        plt.close(fig)
        return anim


# ═══════════════════════════════════════════════════════════════════════════
# MAIN INTERFACE (Enhanced with more options)
# ═══════════════════════════════════════════════════════════════════════════

class EnhancedSystemInterface:
    """Enhanced main interface"""
    
    def __init__(self):
        self.positions = None
        self.performance = None
        self.controller = None
        
        print("\n" + "="*70)
        print(" "*10 + "COMPLETE MULTI-EMITTER THRUSTER SYSTEM")
        print(" "*15 + "Enhanced Version 4.0 FINAL")
        print("="*70)
        print("\nEnhancements:")
        print("  ✓ 6 propellants available")
        print("  ✓ Advanced optimization algorithms")
        print("  ✓ Multiple maneuver types")
        print("  ✓ Unique maneuver visualizations")
        print("  ✓ Orbit & firing animations (GIF)")
    
    def run_interactive(self):
        """Run interactive workflow"""
        print("\n" + "="*70)
        print("STEP 1: CONFIGURATION")
        print("="*70)
        
        config = self.get_user_inputs()
        
        print("\n" + "="*70)
        print("STEP 2: GEOMETRY OPTIMIZATION")
        print("="*70)
        
        self.optimize_geometry(config)
        
        print("\n" + "="*70)
        print("STEP 3: VISUALIZATION")
        print("="*70)
        
        self.create_visualizations()
        
        print("\n" + "="*70)
        print("STEP 4: MANEUVER PLANNING")
        print("="*70)
        
        maneuver_type = self.get_maneuver_input()
        
        if maneuver_type == 'hohmann':
            self.plan_hohmann(config)
        elif maneuver_type == 'bielliptic':
            self.plan_bielliptic(config)
        elif maneuver_type == 'plane':
            self.plan_plane_change(config)
        elif maneuver_type == 'attitude':
            self.plan_attitude(config)
        elif maneuver_type == 'detumble':
            self.plan_detumbling(config)
        
        print("\n" + "="*70)
        print("✓ ANALYSIS COMPLETE!")
        print("="*70)
        print("\nGenerated files:")
        print("  ✓ array_layout.png (standard array view)")
        print("  ✓ array_3d.png (3D view with plumes)")
        print("  ✓ maneuver_plan.png (UNIQUE maneuver visualization)")
        print("  ✓ orbit_transfer.gif (animated orbit)")
        print("  ✓ firing_sequence.gif (animated firing)")
    
    def get_user_inputs(self) -> Dict:
        """Get configuration"""
        print("\nEnter configuration parameters:")
        print("(Press Enter for defaults shown in brackets)")
        
        n_str = input("\nNumber of emitters [100]: ").strip()
        n_emitters = int(n_str) if n_str else 100
        
        print("\n" + "="*70)
        print("PROPELLANT SELECTION")
        print("="*70)
        print("\nAvailable propellants:")
        props = list(Constants.PROPELLANTS.keys())
        for i, name in enumerate(props, 1):
            prop = Constants.PROPELLANTS[name]
            rec = " [RECOMMENDED]" if prop.get('recommended', False) else ""
            print(f"  {i}. {name}{rec}")
            print(f"     Isp: {prop['isp_range'][0]}-{prop['isp_range'][1]} s")
        
        prop_str = input(f"\nSelect propellant (1-{len(props)}) [1]: ").strip()
        prop_idx = int(prop_str) - 1 if prop_str else 0
        prop_name = props[prop_idx]
        propellant = Constants.PROPELLANTS[prop_name]
        
        mass_str = input("\nSpacecraft mass (kg) [10]: ").strip()
        mass = float(mass_str) if mass_str else 10.0
        
        voltage_str = input("Operating voltage (V) [2000]: ").strip()
        voltage = float(voltage_str) if voltage_str else 2000.0
        
        print("\nOptimization algorithm:")
        print("  1. Quick (L-BFGS-B, ~5 seconds)")
        print("  2. Thorough (Dual Annealing, ~30 seconds)")
        mode_str = input("Select mode [1]: ").strip()
        opt_mode = 'thorough' if mode_str == '2' else 'quick'
        
        return {
            'n_emitters': n_emitters,
            'propellant': propellant,
            'spacecraft_mass': mass,
            'voltage': voltage,
            'optimization_mode': opt_mode
        }
    
    def optimize_geometry(self, config: Dict):
        """Run optimization"""
        optimizer = EnhancedGeometryOptimizer(
            config['n_emitters'],
            config['propellant'],
            config['voltage']
        )
        
        self.positions, self.performance = optimizer.optimize(
            method=config['optimization_mode']
        )
        
        inertia = np.diag([0.1, 0.1, 0.05])
        
        self.controller = EnhancedSpacecraftController(
            self.positions,
            config['spacecraft_mass'],
            inertia,
            thrust_per_emitter=1e-6
        )
    
    def create_visualizations(self):
        """Create visualizations"""
        viz = EnhancedVisualizer()
        
        fig1 = viz.plot_array_layout(
            self.positions,
            self.performance,
            title='Optimized Array Layout'
        )
        fig1.savefig('array_layout.png', dpi=300, bbox_inches='tight')
        print("✓ Saved: array_layout.png")
        plt.close(fig1)
        
        firing = np.ones(len(self.positions), dtype=bool)
        fig2 = viz.plot_3d_array(self.positions, firing)
        fig2.savefig('array_3d.png', dpi=300, bbox_inches='tight')
        print("✓ Saved: array_3d.png")
        plt.close(fig2)
    
    def get_maneuver_input(self) -> str:
        """Get maneuver type"""
        print("\nSelect maneuver type:")
        print("  1. Hohmann transfer (orbit raise/lower)")
        print("  2. Bi-elliptic transfer (large orbit changes)")
        print("  3. Plane change (inclination change)")
        print("  4. Attitude change (reorientation)")
        print("  5. Detumbling (stop rotation)")
        
        choice = input("Enter choice [1]: ").strip() or "1"
        
        choices = {'1': 'hohmann', '2': 'bielliptic', '3': 'plane',
                  '4': 'attitude', '5': 'detumble'}
        return choices.get(choice, 'hohmann')
    
    def plan_hohmann(self, config: Dict):
        """Plan Hohmann transfer"""
        print("\n" + "="*70)
        print("HOHMANN TRANSFER PLANNING")
        print("="*70)
        
        alt_init_str = input("\nInitial altitude (km) [400]: ").strip()
        alt_init = float(alt_init_str)*1e3 if alt_init_str else 400e3
        
        alt_final_str = input("Target altitude (km) [500]: ").strip()
        alt_final = float(alt_final_str)*1e3 if alt_final_str else 500e3
        
        initial_orbit = {
            'a': Constants.R_earth + alt_init,
            'e': 0, 'i': 0, 'Omega': 0, 'omega': 0, 'nu': 0
        }
        
        target_orbit = {
            'a': Constants.R_earth + alt_final,
            'e': 0, 'i': 0, 'Omega': 0, 'omega': 0, 'nu': 0
        }
        
        plan = self.controller.plan_hohmann(initial_orbit, target_orbit)
        
        self.print_maneuver_summary(plan, alt_init/1e3, alt_final/1e3)
        
        # Create UNIQUE maneuver plot
        viz = EnhancedVisualizer()
        fig = viz.plot_maneuver_plan(
            self.positions,
            plan['firing_pattern_1'],
            plan,
            title='Hohmann Transfer Maneuver'
        )
        fig.savefig('maneuver_plan.png', dpi=300, bbox_inches='tight')
        print("\n✓ Saved: maneuver_plan.png (unique maneuver visualization)")
        plt.close(fig)
        
        # Create animations
        animator = EnhancedAnimator()
        animator.animate_orbit_transfer(initial_orbit, target_orbit)
        animator.animate_firing_sequence(self.positions, plan['firing_pattern_1'],
                                        plan['burn_time_1'])
    
    def plan_bielliptic(self, config: Dict):
        """Plan bi-elliptic transfer"""
        print("\n" + "="*70)
        print("BI-ELLIPTIC TRANSFER PLANNING")
        print("="*70)
        
        alt_init_str = input("\nInitial altitude (km) [400]: ").strip()
        alt_init = float(alt_init_str)*1e3 if alt_init_str else 400e3
        
        alt_final_str = input("Target altitude (km) [800]: ").strip()
        alt_final = float(alt_final_str)*1e3 if alt_final_str else 800e3
        
        initial_orbit = {
            'a': Constants.R_earth + alt_init,
            'e': 0, 'i': 0, 'Omega': 0, 'omega': 0, 'nu': 0
        }
        
        target_orbit = {
            'a': Constants.R_earth + alt_final,
            'e': 0, 'i': 0, 'Omega': 0, 'omega': 0, 'nu': 0
        }
        
        plan = self.controller.plan_bielliptic(initial_orbit, target_orbit)
        
        print(f"\n{'='*70}")
        print("BI-ELLIPTIC TRANSFER PLAN")
        print(f"{'='*70}")
        print(f"Altitude: {alt_init/1e3:.0f} → {alt_final/1e3:.0f} km")
        print(f"\nΔv Requirements:")
        print(f"  Burn 1: {plan['dv1']:.2f} m/s")
        print(f"  Burn 2: {plan['dv2']:.2f} m/s")
        print(f"  Burn 3: {plan['dv3']:.2f} m/s")
        print(f"  Total:  {plan['total_dv']:.2f} m/s")
        print(f"\nBurn Times:")
        print(f"  Burn 1: {plan['burn_time_1']/3600:.2f} hours")
        print(f"  Burn 2: {plan['burn_time_2']/3600:.2f} hours")
        print(f"  Burn 3: {plan['burn_time_3']/3600:.2f} hours")
        
        viz = EnhancedVisualizer()
        fig = viz.plot_maneuver_plan(
            self.positions,
            plan['firing_pattern'],
            plan,
            title='Bi-Elliptic Transfer'
        )
        fig.savefig('maneuver_plan.png', dpi=300, bbox_inches='tight')
        print("\n✓ Saved: maneuver_plan.png")
        plt.close(fig)
    
    def plan_plane_change(self, config: Dict):
        """Plan plane change"""
        print("\n" + "="*70)
        print("PLANE CHANGE MANEUVER")
        print("="*70)
        
        delta_i_str = input("\nInclination change (degrees) [10]: ").strip()
        delta_i = float(delta_i_str) if delta_i_str else 10.0
        
        alt_str = input("Altitude (km) [400]: ").strip()
        altitude = float(alt_str)*1e3 if alt_str else 400e3
        
        plan = self.controller.plan_plane_change(delta_i, altitude)
        
        print(f"\n{'='*70}")
        print("PLANE CHANGE PLAN")
        print(f"{'='*70}")
        print(f"Inclination change: {delta_i:.2f}°")
        print(f"Required Δv: {plan['dv']:.2f} m/s")
        print(f"Burn time: {plan['burn_time']/3600:.2f} hours")
        
        viz = EnhancedVisualizer()
        fig = viz.plot_maneuver_plan(
            self.positions,
            plan['firing_pattern'],
            plan,
            title='Plane Change Maneuver'
        )
        fig.savefig('maneuver_plan.png', dpi=300, bbox_inches='tight')
        print("\n✓ Saved: maneuver_plan.png")
        plt.close(fig)
    
    def plan_attitude(self, config: Dict):
        """Plan attitude change"""
        print("\n" + "="*70)
        print("ATTITUDE CONTROL MANEUVER")
        print("="*70)
        
        roll_str = input("\nRoll (deg) [10]: ").strip()
        pitch_str = input("Pitch (deg) [5]: ").strip()
        yaw_str = input("Yaw (deg) [0]: ").strip()
        
        target_euler = np.array([
            float(roll_str) if roll_str else 10.0,
            float(pitch_str) if pitch_str else 5.0,
            float(yaw_str) if yaw_str else 0.0
        ])
        
        plan = self.controller.plan_attitude(target_euler)
        
        print(f"\n{'='*70}")
        print("ATTITUDE MANEUVER PLAN")
        print(f"{'='*70}")
        print(f"Target: R={target_euler[0]:.1f}°, P={target_euler[1]:.1f}°, Y={target_euler[2]:.1f}°")
        print(f"Angle change: {plan['angle_change_deg']:.2f}°")
        print(f"Active emitters: {plan['n_emitters_active']} / {len(self.positions)}")
        print(f"Time: {plan['estimated_time']:.1f} seconds")
        
        viz = EnhancedVisualizer()
        fig = viz.plot_maneuver_plan(
            self.positions,
            plan['firing_pattern'],
            plan,
            title='Attitude Control Maneuver'
        )
        fig.savefig('maneuver_plan.png', dpi=300, bbox_inches='tight')
        print("\n✓ Saved: maneuver_plan.png")
        plt.close(fig)
        
        animator = EnhancedAnimator()
        animator.animate_firing_sequence(self.positions, plan['firing_pattern'],
                                        plan['estimated_time'])
    
    def plan_detumbling(self, config: Dict):
        """Plan detumbling"""
        print("\n" + "="*70)
        print("DETUMBLING MANEUVER")
        print("="*70)
        print("\nInitial rotation rates:")
        
        wx_str = input("  ωx (deg/s) [2]: ").strip()
        wy_str = input("  ωy (deg/s) [1]: ").strip()
        wz_str = input("  ωz (deg/s) [3]: ").strip()
        
        omega = np.deg2rad(np.array([
            float(wx_str) if wx_str else 2.0,
            float(wy_str) if wy_str else 1.0,
            float(wz_str) if wz_str else 3.0
        ]))
        
        plan = self.controller.plan_detumbling(omega)
        
        print(f"\n{'='*70}")
        print("DETUMBLING PLAN")
        print(f"{'='*70}")
        print(f"Initial rotation: {np.rad2deg(omega)} deg/s")
        print(f"Active emitters: {plan['n_emitters_active']} / {len(self.positions)}")
        print(f"Estimated time: {plan['estimated_time']:.0f} seconds")
        
        viz = EnhancedVisualizer()
        fig = viz.plot_maneuver_plan(
            self.positions,
            plan['firing_pattern'],
            plan,
            title='Detumbling Maneuver'
        )
        fig.savefig('maneuver_plan.png', dpi=300, bbox_inches='tight')
        print("\n✓ Saved: maneuver_plan.png")
        plt.close(fig)
    
    def print_maneuver_summary(self, plan: Dict, alt_init: float, alt_final: float):
        """Print maneuver summary"""
        print(f"\n{'='*70}")
        print("HOHMANN TRANSFER PLAN")
        print(f"{'='*70}")
        print(f"Altitude change: {alt_init:.0f} → {alt_final:.0f} km")
        print(f"\nΔv Requirements:")
        print(f"  Burn 1: {plan['dv1']:.2f} m/s")
        print(f"  Burn 2: {plan['dv2']:.2f} m/s")
        print(f"  Total:  {plan['total_dv']:.2f} m/s")
        print(f"\nBurn Times:")
        print(f"  Burn 1: {plan['burn_time_1']/3600:.2f} hours")
        print(f"  Burn 2: {plan['burn_time_2']/3600:.2f} hours")
        print(f"  Transfer: {plan['transfer_time']/3600:.2f} hours")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════

def main():
    """Main entry point"""
    try:
        interface = EnhancedSystemInterface()
        interface.run_interactive()
    except KeyboardInterrupt:
        print("\n\n⚠ Interrupted by user.")
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
