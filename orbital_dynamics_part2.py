"""
PART : Complete 6-DOF Spacecraft Dynamics with Thruster Control Allocation
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize, fsolve
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ============================================================================
# ASTRODYNAMICS CONSTANTS
# ============================================================================

class OrbitalConstants:
    """Earth orbital mechanics constants"""
    mu_earth = 3.986004418e14  # Earth gravitational parameter (m³/s²)
    R_earth = 6378137  # Earth equatorial radius (m)
    J2 = 1.08262668e-3  # Earth J2 coefficient (oblateness)
    omega_earth = 7.2921150e-5  # Earth rotation rate (rad/s)
    
    # Standard orbits
    LEO_altitude = 400e3  # Low Earth Orbit (m)
    GEO_altitude = 35786e3  # Geostationary (m)


# ============================================================================
# ORBITAL STATE & ELEMENTS
# ============================================================================

class OrbitalState:
    """
    Complete spacecraft orbital state
    
    Can represent state as:
    - Cartesian (position, velocity in ECI)
    - Classical Orbital Elements (COE)
    - Equinoctial Elements
    """
    
    def __init__(self, r: np.ndarray, v: np.ndarray, t: float = 0):
        """
        Args:
            r: Position in ECI (m) [x, y, z]
            v: Velocity in ECI (m/s) [vx, vy, vz]
            t: Time (s)
        """
        self.r = np.array(r, dtype=float)
        self.v = np.array(v, dtype=float)
        self.t = t
        
        # Compute derived quantities
        self._compute_derived()
    
    def _compute_derived(self):
        """Compute orbital elements and other derived quantities"""
        r_mag = np.linalg.norm(self.r)
        v_mag = np.linalg.norm(self.v)
        
        # Specific orbital energy
        self.energy = v_mag**2 / 2 - OrbitalConstants.mu_earth / r_mag
        
        # Angular momentum vector
        self.h = np.cross(self.r, self.v)
        h_mag = np.linalg.norm(self.h)
        
        # Eccentricity vector
        self.e_vec = ((v_mag**2 - OrbitalConstants.mu_earth/r_mag) * self.r - 
                     np.dot(self.r, self.v) * self.v) / OrbitalConstants.mu_earth
        self.e = np.linalg.norm(self.e_vec)  # Eccentricity
        
        # Semi-major axis
        if abs(self.e - 1.0) > 1e-10:  # Not parabolic
            self.a = -OrbitalConstants.mu_earth / (2 * self.energy)
        else:
            self.a = np.inf
        
        # Inclination
        self.i = np.arccos(self.h[2] / h_mag)
        
        # Node vector
        N = np.cross([0, 0, 1], self.h)
        N_mag = np.linalg.norm(N)
        
        # Right ascension of ascending node (RAAN)
        if N_mag > 1e-10:
            self.Omega = np.arctan2(N[1], N[0])
            if self.Omega < 0:
                self.Omega += 2 * np.pi
        else:
            self.Omega = 0
        
        # Argument of periapsis
        if N_mag > 1e-10 and self.e > 1e-10:
            cos_omega = np.dot(N, self.e_vec) / (N_mag * self.e)
            cos_omega = np.clip(cos_omega, -1, 1)
            self.omega = np.arccos(cos_omega)
            if self.e_vec[2] < 0:
                self.omega = 2*np.pi - self.omega
        else:
            self.omega = 0
        
        # True anomaly
        if self.e > 1e-10:
            cos_nu = np.dot(self.e_vec, self.r) / (self.e * r_mag)
            cos_nu = np.clip(cos_nu, -1, 1)
            self.nu = np.arccos(cos_nu)
            if np.dot(self.r, self.v) < 0:
                self.nu = 2*np.pi - self.nu
        else:
            self.nu = 0
        
        # Period (for elliptical orbits)
        if self.a > 0:
            self.T = 2 * np.pi * np.sqrt(self.a**3 / OrbitalConstants.mu_earth)
        else:
            self.T = np.inf
    
    def get_classical_elements(self) -> Dict:
        """Return classical orbital elements"""
        return {
            'a': self.a,  # Semi-major axis (m)
            'e': self.e,  # Eccentricity
            'i': self.i,  # Inclination (rad)
            'Omega': self.Omega,  # RAAN (rad)
            'omega': self.omega,  # Argument of periapsis (rad)
            'nu': self.nu,  # True anomaly (rad)
            'T': self.T  # Orbital period (s)
        }
    
    def altitude(self) -> float:
        """Altitude above Earth surface (m)"""
        return np.linalg.norm(self.r) - OrbitalConstants.R_earth
    
    @classmethod
    def from_coe(cls, a: float, e: float, i: float, Omega: float,
                omega: float, nu: float):
        """
        Create from classical orbital elements
        
        Args:
            a: Semi-major axis (m)
            e: Eccentricity
            i: Inclination (rad)
            Omega: RAAN (rad)
            omega: Argument of periapsis (rad)
            nu: True anomaly (rad)
        """
        # Orbital parameter
        p = a * (1 - e**2)
        
        # Position and velocity in perifocal frame
        r_pf = p / (1 + e*np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0])
        v_pf = np.sqrt(OrbitalConstants.mu_earth/p) * np.array([
            -np.sin(nu),
            e + np.cos(nu),
            0
        ])
        
        # Rotation matrix perifocal -> ECI
        R3_Omega = np.array([
            [np.cos(Omega), -np.sin(Omega), 0],
            [np.sin(Omega), np.cos(Omega), 0],
            [0, 0, 1]
        ])
        
        R1_i = np.array([
            [1, 0, 0],
            [0, np.cos(i), -np.sin(i)],
            [0, np.sin(i), np.cos(i)]
        ])
        
        R3_omega = np.array([
            [np.cos(omega), -np.sin(omega), 0],
            [np.sin(omega), np.cos(omega), 0],
            [0, 0, 1]
        ])
        
        R_pf2eci = R3_Omega @ R1_i @ R3_omega
        
        # Transform to ECI
        r = R_pf2eci @ r_pf
        v = R_pf2eci @ v_pf
        
        return cls(r, v)


# ============================================================================
# ORBITAL PROPAGATION
# ============================================================================

class OrbitalPropagator:
    """
    Propagate spacecraft orbit under various force models
    
    Force models:
    - Two-body (Keplerian)
    - J2 perturbation (Earth oblateness)
    - Atmospheric drag
    - Solar radiation pressure
    - Third-body (Moon, Sun)
    - Thrust
    """
    
    def __init__(self, force_models: List[str] = ['two_body']):
        """
        Args:
            force_models: List of force models to include
        """
        self.force_models = force_models
        
        print(f"Orbital Propagator initialized with models: {force_models}")
    
    def equations_of_motion(self, t: float, state: np.ndarray,
                           thrust_func=None, spacecraft_params: Dict = None) -> np.ndarray:
        """
        Equations of motion: ṙ = v, v̇ = a
        
        Args:
            t: Time (s)
            state: [x, y, z, vx, vy, vz]
            thrust_func: Function returning thrust vector at time t
            spacecraft_params: Mass, area, etc.
        
        Returns:
            state_dot: [vx, vy, vz, ax, ay, az]
        """
        r = state[0:3]
        v = state[3:6]
        
        r_mag = np.linalg.norm(r)
        
        # Acceleration
        a = np.zeros(3)
        
        # Two-body (Keplerian)
        if 'two_body' in self.force_models:
            a += -OrbitalConstants.mu_earth * r / r_mag**3
        
        # J2 perturbation
        if 'J2' in self.force_models:
            z2 = r[2]**2
            r2 = r_mag**2
            R_E = OrbitalConstants.R_earth
            J2 = OrbitalConstants.J2
            
            factor = 3/2 * J2 * OrbitalConstants.mu_earth * R_E**2 / r_mag**5
            
            a[0] += factor * r[0] * (5*z2/r2 - 1)
            a[1] += factor * r[1] * (5*z2/r2 - 1)
            a[2] += factor * r[2] * (5*z2/r2 - 3)
        
        # Atmospheric drag
        if 'drag' in self.force_models and spacecraft_params:
            # Simplified exponential atmosphere
            h = r_mag - OrbitalConstants.R_earth
            rho0 = 1.225  # kg/m³ at sea level
            H = 8500  # Scale height (m)
            rho = rho0 * np.exp(-h / H)
            
            # Drag acceleration
            Cd = spacecraft_params.get('Cd', 2.2)
            A = spacecraft_params.get('area', 0.01)  # m²
            m = spacecraft_params.get('mass', 10)  # kg
            
            v_rel = v  # Relative to atmosphere (simplified)
            v_rel_mag = np.linalg.norm(v_rel)
            
            if v_rel_mag > 0:
                a_drag = -0.5 * rho * Cd * A / m * v_rel_mag * v_rel
                a += a_drag
        
        # Thrust
        if thrust_func is not None:
            thrust_vec = thrust_func(t)  # Thrust in inertial frame
            m = spacecraft_params.get('mass', 10) if spacecraft_params else 10
            a += thrust_vec / m
        
        # State derivative
        state_dot = np.concatenate([v, a])
        
        return state_dot
    
    def propagate(self, initial_state: OrbitalState, duration: float,
                 thrust_func=None, spacecraft_params: Dict = None,
                 dt: float = 60) -> Tuple[np.ndarray, np.ndarray]:
        """
        Propagate orbit
        
        Args:
            initial_state: Initial orbital state
            duration: Propagation time (s)
            thrust_func: Optional thrust function
            spacecraft_params: Spacecraft parameters
            dt: Output time step (s)
        
        Returns:
            (t_array, state_array): Time and state history
        """
        # Initial state vector
        y0 = np.concatenate([initial_state.r, initial_state.v])
        
        # Time span
        t_span = (0, duration)
        t_eval = np.arange(0, duration, dt)
        
        # Integrate
        sol = solve_ivp(
            lambda t, y: self.equations_of_motion(t, y, thrust_func, spacecraft_params),
            t_span,
            y0,
            t_eval=t_eval,
            method='DOP853',  # High-order Runge-Kutta
            rtol=1e-9,
            atol=1e-12
        )
        
        return sol.t, sol.y.T


# ============================================================================
# ATTITUDE DYNAMICS
# ============================================================================

class AttitudeDynamics:
    """
    Spacecraft attitude dynamics using quaternions
    
    Euler's equations + quaternion kinematics
    """
    
    def __init__(self, inertia_matrix: np.ndarray):
        """
        Args:
            inertia_matrix: 3x3 moment of inertia matrix (kg·m²)
        """
        self.I = inertia_matrix
        self.I_inv = np.linalg.inv(inertia_matrix)
        
        print(f"Attitude Dynamics initialized:")
        print(f"  Inertia matrix (kg·m²):")
        print(f"    {self.I}")
    
    def quaternion_multiply(self, q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
        """
        Quaternion multiplication: q1 ⊗ q2
        
        Quaternion format: [q0, q1, q2, q3] where q0 is scalar part
        """
        w1, x1, y1, z1 = q1
        w2, x2, y2, z2 = q2
        
        return np.array([
            w1*w2 - x1*x2 - y1*y2 - z1*z2,
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2
        ])
    
    def quaternion_derivative(self, q: np.ndarray, omega: np.ndarray) -> np.ndarray:
        """
        Quaternion kinematic equation: q̇ = 0.5 * Ω(ω) * q
        
        Args:
            q: Quaternion [q0, q1, q2, q3]
            omega: Angular velocity in body frame (rad/s)
        
        Returns:
            q_dot: Quaternion derivative
        """
        wx, wy, wz = omega
        
        # Omega matrix
        Omega = np.array([
            [0, -wx, -wy, -wz],
            [wx, 0, wz, -wy],
            [wy, -wz, 0, wx],
            [wz, wy, -wx, 0]
        ])
        
        q_dot = 0.5 * Omega @ q
        
        return q_dot
    
    def eulers_equation(self, omega: np.ndarray, torque: np.ndarray) -> np.ndarray:
        """
        Euler's rotational equations of motion
        
        I·ω̇ + ω × (I·ω) = τ
        
        Args:
            omega: Angular velocity (rad/s)
            torque: Applied torque (N·m)
        
        Returns:
            omega_dot: Angular acceleration (rad/s²)
        """
        I_omega = self.I @ omega
        omega_cross_I_omega = np.cross(omega, I_omega)
        
        omega_dot = self.I_inv @ (torque - omega_cross_I_omega)
        
        return omega_dot
    
    def propagate_attitude(self, q0: np.ndarray, omega0: np.ndarray,
                          duration: float, torque_func,
                          dt: float = 0.1) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Propagate attitude dynamics
        
        Args:
            q0: Initial quaternion
            omega0: Initial angular velocity (rad/s)
            duration: Time (s)
            torque_func: Function returning torque at time t
            dt: Time step (s)
        
        Returns:
            (t, q, omega): Time, quaternion history, angular velocity history
        """
        def state_derivative(t, state):
            q = state[0:4]
            omega = state[4:7]
            
            # Normalize quaternion
            q = q / np.linalg.norm(q)
            
            # Get torque
            torque = torque_func(t)
            
            # Derivatives
            q_dot = self.quaternion_derivative(q, omega)
            omega_dot = self.eulers_equation(omega, torque)
            
            return np.concatenate([q_dot, omega_dot])
        
        # Initial state
        y0 = np.concatenate([q0, omega0])
        
        # Time span
        t_span = (0, duration)
        t_eval = np.arange(0, duration, dt)
        
        # Integrate
        sol = solve_ivp(
            state_derivative,
            t_span,
            y0,
            t_eval=t_eval,
            method='RK45'
        )
        
        # Extract results
        q_history = sol.y[0:4, :].T
        omega_history = sol.y[4:7, :].T
        
        # Normalize quaternions
        for i in range(len(q_history)):
            q_history[i] = q_history[i] / np.linalg.norm(q_history[i])
        
        return sol.t, q_history, omega_history
    
    @staticmethod
    def quaternion_to_dcm(q: np.ndarray) -> np.ndarray:
        """
        Convert quaternion to Direction Cosine Matrix (DCM)
        
        DCM rotates vectors from body frame to inertial frame
        """
        q0, q1, q2, q3 = q
        
        DCM = np.array([
            [q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]
        ])
        
        return DCM
    
    @staticmethod
    def quaternion_to_euler(q: np.ndarray) -> np.ndarray:
        """
        Convert quaternion to Euler angles (3-2-1 sequence)
        
        Returns: [roll, pitch, yaw] in radians
        """
        q0, q1, q2, q3 = q
        
        # Roll (φ)
        roll = np.arctan2(2*(q0*q1 + q2*q3), 1 - 2*(q1**2 + q2**2))
        
        # Pitch (θ)
        sin_pitch = 2*(q0*q2 - q3*q1)
        sin_pitch = np.clip(sin_pitch, -1, 1)
        pitch = np.arcsin(sin_pitch)
        
        # Yaw (ψ)
        yaw = np.arctan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
        
        return np.array([roll, pitch, yaw])


# ============================================================================
# TO BE CONTINUED IN PART 3 (CONTROL ALLOCATION)
# ============================================================================
