"""
Configuration and Validation System for Electrospray Thruster Simulation
==========================================================================

Provides:
- Parameter validation with physical bounds
- Configuration management
- Unit conversions
- Physical constants
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Any
import warnings
import json


# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

class PhysicalConstants:
    """Universal physical constants"""
    EPSILON_0 = 8.854187817e-12  # F/m - Permittivity of free space
    E_CHARGE = 1.602176634e-19   # C - Elementary charge
    K_B = 1.380649e-23           # J/K - Boltzmann constant
    G_0 = 9.81                   # m/s² - Standard gravity
    AMU = 1.66053906660e-27      # kg - Atomic mass unit
    N_A = 6.02214076e23          # Avogadro's number
    R_GAS = 8.314                # J/(mol·K) - Universal gas constant
    
    @staticmethod
    def get_all() -> Dict[str, float]:
        """Return all constants as dictionary"""
        return {
            'EPSILON_0': PhysicalConstants.EPSILON_0,
            'E_CHARGE': PhysicalConstants.E_CHARGE,
            'K_B': PhysicalConstants.K_B,
            'G_0': PhysicalConstants.G_0,
            'AMU': PhysicalConstants.AMU,
            'N_A': PhysicalConstants.N_A,
            'R_GAS': PhysicalConstants.R_GAS
        }


# ============================================================================
# UNIT CONVERSIONS
# ============================================================================

class UnitConverter:
    """Handle unit conversions"""
    
    @staticmethod
    def voltage(value: float, from_unit: str = 'V', to_unit: str = 'V') -> float:
        """Convert voltage units"""
        conv = {
            'V': 1.0,
            'kV': 1e3,
            'mV': 1e-3,
            'μV': 1e-6
        }
        return value * conv[from_unit] / conv[to_unit]
    
    @staticmethod
    def length(value: float, from_unit: str = 'm', to_unit: str = 'm') -> float:
        """Convert length units"""
        conv = {
            'm': 1.0,
            'mm': 1e-3,
            'μm': 1e-6,
            'nm': 1e-9,
            'Å': 1e-10
        }
        return value * conv[from_unit] / conv[to_unit]
    
    @staticmethod
    def mass_flow(value: float, from_unit: str = 'kg/s', to_unit: str = 'kg/s') -> float:
        """Convert mass flow rate units"""
        conv = {
            'kg/s': 1.0,
            'g/s': 1e-3,
            'mg/s': 1e-6,
            'μg/s': 1e-9,
            'ng/s': 1e-12,
            'pg/s': 1e-15
        }
        return value * conv[from_unit] / conv[to_unit]
    
    @staticmethod
    def current(value: float, from_unit: str = 'A', to_unit: str = 'A') -> float:
        """Convert current units"""
        conv = {
            'A': 1.0,
            'mA': 1e-3,
            'μA': 1e-6,
            'nA': 1e-9,
            'pA': 1e-12
        }
        return value * conv[from_unit] / conv[to_unit]
    
    @staticmethod
    def thrust(value: float, from_unit: str = 'N', to_unit: str = 'N') -> float:
        """Convert thrust units"""
        conv = {
            'N': 1.0,
            'mN': 1e-3,
            'μN': 1e-6,
            'nN': 1e-9
        }
        return value * conv[from_unit] / conv[to_unit]


# ============================================================================
# PARAMETER VALIDATION
# ============================================================================

@dataclass
class ParameterBounds:
    """Physical bounds for simulation parameters"""
    
    # Voltage bounds (Volts)
    V_min: float = 100.0
    V_max: float = 20000.0
    
    # Mass flow rate bounds (kg/s)
    m_dot_min: float = 1e-15
    m_dot_max: float = 1e-6
    
    # Geometric bounds (meters)
    d_tip_min: float = 1e-6
    d_tip_max: float = 1e-3
    gap_min: float = 10e-6
    gap_max: float = 50e-3
    
    # Temperature bounds (Kelvin)
    T_min: float = 200.0
    T_max: float = 600.0
    
    # Liquid property bounds
    conductivity_min: float = 0.01  # S/m
    conductivity_max: float = 100.0
    surface_tension_min: float = 0.01  # N/m
    surface_tension_max: float = 0.2
    viscosity_min: float = 1e-4  # Pa·s
    viscosity_max: float = 1.0
    density_min: float = 500.0  # kg/m³
    density_max: float = 3000.0
    
    def validate_voltage(self, V: float, name: str = "Voltage") -> None:
        """Validate voltage parameter"""
        if not (self.V_min <= V <= self.V_max):
            raise ValueError(
                f"{name} = {V:.1f} V is outside valid range "
                f"[{self.V_min:.1f}, {self.V_max:.1f}] V"
            )
    
    def validate_mass_flow(self, m_dot: float, name: str = "Mass flow") -> None:
        """Validate mass flow rate"""
        if not (self.m_dot_min <= m_dot <= self.m_dot_max):
            raise ValueError(
                f"{name} = {m_dot:.2e} kg/s is outside valid range "
                f"[{self.m_dot_min:.2e}, {self.m_dot_max:.2e}] kg/s"
            )
    
    def validate_geometry(self, d: float, name: str = "Dimension") -> None:
        """Validate geometric parameter"""
        if not (self.d_tip_min <= d <= self.d_tip_max):
            raise ValueError(
                f"{name} = {d:.2e} m is outside valid range "
                f"[{self.d_tip_min:.2e}, {self.d_tip_max:.2e}] m"
            )
    
    def validate_temperature(self, T: float, name: str = "Temperature") -> None:
        """Validate temperature"""
        if not (self.T_min <= T <= self.T_max):
            raise ValueError(
                f"{name} = {T:.1f} K is outside valid range "
                f"[{self.T_min:.1f}, {self.T_max:.1f}] K"
            )
    
    def validate_liquid_property(self, value: float, prop_type: str) -> None:
        """Validate liquid property"""
        bounds_map = {
            'conductivity': (self.conductivity_min, self.conductivity_max, 'S/m'),
            'surface_tension': (self.surface_tension_min, self.surface_tension_max, 'N/m'),
            'viscosity': (self.viscosity_min, self.viscosity_max, 'Pa·s'),
            'density': (self.density_min, self.density_max, 'kg/m³')
        }
        
        if prop_type not in bounds_map:
            raise ValueError(f"Unknown liquid property type: {prop_type}")
        
        min_val, max_val, unit = bounds_map[prop_type]
        if not (min_val <= value <= max_val):
            raise ValueError(
                f"{prop_type} = {value:.3e} {unit} is outside valid range "
                f"[{min_val:.3e}, {max_val:.3e}] {unit}"
            )


# ============================================================================
# SIMULATION CONFIGURATION
# ============================================================================

@dataclass
class SimulationConfig:
    """Configuration for electrospray simulation"""
    
    # Validation
    validate_inputs: bool = True
    warn_on_unstable: bool = True
    
    # Physics models
    use_self_heating: bool = True
    use_polydispersity: bool = False
    use_space_charge: bool = False
    
    # Numerical parameters
    convergence_tolerance: float = 1e-6
    max_iterations: int = 1000
    
    # Output control
    verbose: bool = True
    save_intermediate: bool = False
    
    # Temperature model
    ambient_temperature: float = 298.15  # K
    
    # Parameter bounds
    bounds: ParameterBounds = field(default_factory=ParameterBounds)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'validate_inputs': self.validate_inputs,
            'warn_on_unstable': self.warn_on_unstable,
            'use_self_heating': self.use_self_heating,
            'use_polydispersity': self.use_polydispersity,
            'use_space_charge': self.use_space_charge,
            'convergence_tolerance': self.convergence_tolerance,
            'max_iterations': self.max_iterations,
            'verbose': self.verbose,
            'save_intermediate': self.save_intermediate,
            'ambient_temperature': self.ambient_temperature
        }
    
    def save(self, filepath: str) -> None:
        """Save configuration to JSON file"""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    @classmethod
    def load(cls, filepath: str) -> 'SimulationConfig':
        """Load configuration from JSON file"""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls(**data)
    
    def validate_config(self) -> Tuple[bool, list]:
        """Validate configuration consistency"""
        issues = []
        
        if self.convergence_tolerance <= 0:
            issues.append("Convergence tolerance must be positive")
        
        if self.max_iterations < 10:
            issues.append("Max iterations should be at least 10")
        
        if self.ambient_temperature < 100 or self.ambient_temperature > 1000:
            issues.append(f"Ambient temperature {self.ambient_temperature} K seems unrealistic")
        
        return (len(issues) == 0), issues


# ============================================================================
# VALIDATION UTILITIES
# ============================================================================

class ValidationError(Exception):
    """Custom exception for validation errors"""
    pass


class WarningManager:
    """Manage warnings for simulation"""
    
    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self.warnings_issued = []
    
    def warn(self, message: str, category: str = "General") -> None:
        """Issue a warning"""
        if self.enabled:
            warning_msg = f"[{category}] {message}"
            warnings.warn(warning_msg, UserWarning)
            self.warnings_issued.append((category, message))
    
    def get_warnings(self) -> list:
        """Get all warnings issued"""
        return self.warnings_issued
    
    def clear(self) -> None:
        """Clear warning history"""
        self.warnings_issued = []


def check_physical_consistency(params: Dict[str, float]) -> Tuple[bool, str]:
    """
    Check physical consistency of parameters
    
    Returns: (is_consistent, message)
    """
    # Check jet stability (Weber number)
    if 'Q' in params and 'R_jet' in params:
        gamma = params.get('surface_tension', 0.05)
        rho = params.get('density', 1000)
        Q = params['Q']
        R = params['R_jet']
        
        We = rho * Q**2 / (2 * np.pi * R**3 * gamma)
        if We > 1:
            return False, f"Weber number {We:.2f} > 1: jet unstable!"
    
    # Check Reynolds number for viscous stability
    if all(k in params for k in ['surface_tension', 'density', 'viscosity', 'conductivity']):
        gamma = params['surface_tension']
        rho = params['density']
        mu = params['viscosity']
        K = params['conductivity']
        
        Re_K = ((gamma**2 * rho * PhysicalConstants.EPSILON_0) / 
                (mu**3 * K))**(1/3)
        
        if Re_K < 1:
            return False, f"EHD Reynolds {Re_K:.2f} < 1: viscous instability likely!"
    
    return True, "Parameters are physically consistent"


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Configuration and Validation System")
    print("=" * 60)
    
    # Create default config
    config = SimulationConfig(
        use_self_heating=True,
        use_polydispersity=False,
        verbose=True
    )
    
    print("\nDefault Configuration:")
    print(json.dumps(config.to_dict(), indent=2))
    
    # Validate config
    is_valid, issues = config.validate_config()
    if is_valid:
        print("\n✓ Configuration is valid")
    else:
        print("\n✗ Configuration has issues:")
        for issue in issues:
            print(f"  - {issue}")
    
    # Test parameter validation
    bounds = ParameterBounds()
    
    print("\nParameter Validation Tests:")
    print("-" * 60)
    
    try:
        bounds.validate_voltage(2000, "Test voltage")
        print("✓ Voltage 2000 V: Valid")
    except ValueError as e:
        print(f"✗ {e}")
    
    try:
        bounds.validate_voltage(50000, "Too high voltage")
        print("✓ Voltage 50000 V: Valid")
    except ValueError as e:
        print(f"✗ {e}")
    
    # Test unit conversions
    print("\nUnit Conversion Tests:")
    print("-" * 60)
    
    V_kV = UnitConverter.voltage(2000, 'V', 'kV')
    print(f"2000 V = {V_kV} kV")
    
    m_flow_pg = UnitConverter.mass_flow(5e-11, 'kg/s', 'pg/s')
    print(f"5e-11 kg/s = {m_flow_pg} pg/s")
    
    d_um = UnitConverter.length(20e-6, 'm', 'μm')
    print(f"20e-6 m = {d_um} μm")
    
    print("\n" + "=" * 60)
