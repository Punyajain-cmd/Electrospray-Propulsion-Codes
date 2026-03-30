"""
Ionic Liquid Property Database
================================

Comprehensive database of ionic liquid properties for electrospray simulation.
Includes:
- Standard room-temperature ionic liquids
- Temperature-dependent properties
- Validation and uncertainty estimates
- Property interpolation
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List
import warnings

from config import PhysicalConstants, ParameterBounds


@dataclass
class IonicLiquid:
    """
    Ionic liquid propellant properties
    
    All properties at 298.15 K unless otherwise specified
    """
    name: str
    
    # Core properties (at reference temperature)
    conductivity: float        # S/m
    surface_tension: float     # N/m
    viscosity: float           # Pa·s
    density: float             # kg/m³
    permittivity: float        # relative
    
    # Molecular weights
    M_cation: float            # kg/mol
    M_anion: float             # kg/mol
    
    # Temperature dependence (Arrhenius-type)
    K_0: Optional[float] = None          # Pre-exponential for conductivity
    E_a_K: Optional[float] = None        # Activation energy for conductivity (J/mol)
    mu_0: Optional[float] = None         # Pre-exponential for viscosity
    E_a_mu: Optional[float] = None       # Activation energy for viscosity (J/mol)
    
    # Property uncertainties (relative, e.g., 0.1 = 10%)
    conductivity_uncertainty: float = 0.1
    surface_tension_uncertainty: float = 0.05
    viscosity_uncertainty: float = 0.15
    density_uncertainty: float = 0.02
    
    # Reference temperature
    T_ref: float = 298.15  # K
    
    # Additional properties
    electrochemical_window: Optional[float] = None  # V
    thermal_stability: Optional[float] = None       # K (max temperature)
    
    # Metadata
    cas_number: Optional[str] = None
    formula_cation: Optional[str] = None
    formula_anion: Optional[str] = None
    references: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        """Initialize temperature coefficients if not provided"""
        if self.K_0 is None or self.E_a_K is None:
            # Typical activation energy ~0.3 eV for ionic liquids
            self.E_a_K = 0.3 * PhysicalConstants.E_CHARGE * PhysicalConstants.N_A  # J/mol
            self.K_0 = self.conductivity * np.exp(self.E_a_K / (PhysicalConstants.R_GAS * self.T_ref))
        
        if self.mu_0 is None or self.E_a_mu is None:
            # Typical activation energy ~0.4 eV for viscosity
            self.E_a_mu = 0.4 * PhysicalConstants.E_CHARGE * PhysicalConstants.N_A  # J/mol
            self.mu_0 = self.viscosity * np.exp(-self.E_a_mu / (PhysicalConstants.R_GAS * self.T_ref))
        
        # Validate properties
        self._validate_properties()
    
    def _validate_properties(self):
        """Validate all properties are within physical bounds"""
        bounds = ParameterBounds()
        
        try:
            bounds.validate_liquid_property(self.conductivity, 'conductivity')
            bounds.validate_liquid_property(self.surface_tension, 'surface_tension')
            bounds.validate_liquid_property(self.viscosity, 'viscosity')
            bounds.validate_liquid_property(self.density, 'density')
        except ValueError as e:
            warnings.warn(f"Property validation warning for {self.name}: {e}")
    
    def conductivity_T(self, T: float) -> float:
        """Temperature-dependent conductivity using Arrhenius equation"""
        return self.K_0 * np.exp(-self.E_a_K / (PhysicalConstants.R_GAS * T))
    
    def viscosity_T(self, T: float) -> float:
        """Temperature-dependent viscosity using Arrhenius equation"""
        return self.mu_0 * np.exp(self.E_a_mu / (PhysicalConstants.R_GAS * T))
    
    def surface_tension_T(self, T: float) -> float:
        """
        Temperature-dependent surface tension
        Uses linear approximation: γ(T) = γ₀(1 - β(T - T₀))
        Typical β ~ 1×10⁻⁴ K⁻¹
        """
        beta = 1e-4  # K⁻¹
        return self.surface_tension * (1 - beta * (T - self.T_ref))
    
    def density_T(self, T: float) -> float:
        """
        Temperature-dependent density
        Linear thermal expansion: ρ(T) = ρ₀(1 - α(T - T₀))
        Typical α ~ 6×10⁻⁴ K⁻¹
        """
        alpha = 6e-4  # K⁻¹
        return self.density * (1 - alpha * (T - self.T_ref))
    
    def get_properties_at_T(self, T: float) -> Dict[str, float]:
        """Get all temperature-dependent properties"""
        return {
            'conductivity': self.conductivity_T(T),
            'viscosity': self.viscosity_T(T),
            'surface_tension': self.surface_tension_T(T),
            'density': self.density_T(T),
            'permittivity': self.permittivity,  # Assumed constant
            'temperature': T
        }
    
    def average_molecular_weight(self) -> float:
        """Average molecular weight of ion pair"""
        return (self.M_cation + self.M_anion) / 2
    
    def __repr__(self) -> str:
        return (f"IonicLiquid(name='{self.name}', "
                f"σ={self.conductivity:.2f} S/m, "
                f"γ={self.surface_tension:.3f} N/m, "
                f"μ={self.viscosity*1e3:.2f} mPa·s, "
                f"ρ={self.density:.0f} kg/m³)")


# ============================================================================
# STANDARD IONIC LIQUIDS DATABASE
# ============================================================================

# EMI-Im (1-Ethyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide)
EMI_IM = IonicLiquid(
    name="EMI-Im",
    conductivity=0.92,                  # S/m
    surface_tension=0.036,              # N/m
    viscosity=0.0396,                   # Pa·s (39.6 mPa·s)
    density=1520.0,                     # kg/m³
    permittivity=13.8,
    M_cation=111.2e-3,                  # kg/mol
    M_anion=280.1e-3,                   # kg/mol
    electrochemical_window=4.5,         # V
    thermal_stability=673,              # K (400°C)
    cas_number="174899-82-2",
    formula_cation="C6H11N2+",
    formula_anion="C2F6NO4S2-",
    references=[
        "Lozano, P. C., J. Phys. D: Appl. Phys. 39, 126 (2006)",
        "Caballero-Pérez, R. P., J. Appl. Phys. 137, 013301 (2025)"
    ]
)

# EMI-BF4 (1-Ethyl-3-methylimidazolium tetrafluoroborate)
EMI_BF4 = IonicLiquid(
    name="EMI-BF4",
    conductivity=1.46,                  # S/m
    surface_tension=0.072,              # N/m
    viscosity=0.0011,                   # Pa·s (1.1 mPa·s) - Low viscosity
    density=998.0,                      # kg/m³
    permittivity=80.0,
    M_cation=111.17e-3,                 # kg/mol
    M_anion=86.81e-3,                   # kg/mol
    electrochemical_window=4.3,         # V
    thermal_stability=573,              # K (300°C)
    cas_number="143314-16-3",
    formula_cation="C6H11N2+",
    formula_anion="BF4-",
    conductivity_uncertainty=0.08,
    references=[
        "Various electrochemical studies",
        "Standard IL database"
    ]
)

# EAN (Ethylammonium nitrate)
EAN = IonicLiquid(
    name="EAN",
    conductivity=2.27,                  # S/m - High conductivity
    surface_tension=0.0481,             # N/m
    viscosity=0.0398,                   # Pa·s (39.8 mPa·s)
    density=1210.0,                     # kg/m³
    permittivity=29.0,
    M_cation=46.1e-3,                   # kg/mol - Light cation
    M_anion=62.0e-3,                    # kg/mol
    electrochemical_window=3.5,         # V
    thermal_stability=523,              # K (250°C)
    cas_number="22113-86-6",
    formula_cation="C2H8N+",
    formula_anion="NO3-",
    references=[
        "Caballero-Pérez, R. P., J. Appl. Phys. 137, 013301 (2025)",
        "Protic ionic liquid studies"
    ]
)

# BMI-TCM (1-Butyl-3-methylimidazolium tricyanomethanide)
BMI_TCM = IonicLiquid(
    name="BMI-TCM",
    conductivity=1.03,                  # S/m
    surface_tension=0.0496,             # N/m
    viscosity=0.0257,                   # Pa·s (25.7 mPa·s)
    density=1050.0,                     # kg/m³
    permittivity=14.0,
    M_cation=139.2e-3,                  # kg/mol
    M_anion=90.1e-3,                    # kg/mol
    electrochemical_window=4.0,         # V
    thermal_stability=623,              # K (350°C)
    formula_cation="C8H15N2+",
    formula_anion="C4N3-",
    references=[
        "Caballero-Pérez, R. P., J. Appl. Phys. 137, 013301 (2025)"
    ]
)

# BMI-IM (1-Butyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide)
BMI_IM = IonicLiquid(
    name="BMI-IM",
    conductivity=0.38,                  # S/m
    surface_tension=0.0373,             # N/m
    viscosity=0.0524,                   # Pa·s (52.4 mPa·s)
    density=1430.0,                     # kg/m³
    permittivity=11.5,
    M_cation=139.2e-3,                  # kg/mol
    M_anion=280.1e-3,                   # kg/mol
    electrochemical_window=4.2,         # V
    thermal_stability=673,              # K (400°C)
    cas_number="174501-64-5",
    formula_cation="C8H15N2+",
    formula_anion="C2F6NO4S2-",
    references=[
        "Standard IL database",
        "Electrochemical measurements"
    ]
)

# EMIM-DCA (1-Ethyl-3-methylimidazolium dicyanamide)
EMIM_DCA = IonicLiquid(
    name="EMIM-DCA",
    conductivity=1.85,                  # S/m - Very high conductivity
    surface_tension=0.0478,             # N/m
    viscosity=0.0216,                   # Pa·s (21.6 mPa·s)
    density=1100.0,                     # kg/m³
    permittivity=15.0,
    M_cation=111.2e-3,                  # kg/mol
    M_anion=66.04e-3,                   # kg/mol - Light anion
    electrochemical_window=3.8,         # V
    thermal_stability=573,              # K (300°C)
    formula_cation="C6H11N2+",
    formula_anion="C2N3-",
    references=[
        "High-performance IL studies"
    ]
)


# ============================================================================
# IONIC LIQUID DATABASE CLASS
# ============================================================================

class IonicLiquidDatabase:
    """Database of ionic liquid properties"""
    
    def __init__(self):
        self.liquids = {
            'EMI-Im': EMI_IM,
            'EMI-BF4': EMI_BF4,
            'EAN': EAN,
            'BMI-TCM': BMI_TCM,
            'BMI-IM': BMI_IM,
            'EMIM-DCA': EMIM_DCA
        }
    
    def get(self, name: str) -> IonicLiquid:
        """Get liquid by name"""
        if name not in self.liquids:
            available = ', '.join(self.liquids.keys())
            raise ValueError(
                f"Liquid '{name}' not in database. "
                f"Available: {available}"
            )
        return self.liquids[name]
    
    def list_all(self) -> List[str]:
        """List all available liquids"""
        return list(self.liquids.keys())
    
    def add_liquid(self, liquid: IonicLiquid) -> None:
        """Add a new liquid to the database"""
        if liquid.name in self.liquids:
            warnings.warn(
                f"Liquid '{liquid.name}' already exists. Overwriting."
            )
        self.liquids[liquid.name] = liquid
    
    def compare_liquids(self, names: List[str] = None) -> Dict:
        """Compare properties of multiple liquids"""
        if names is None:
            names = self.list_all()
        
        comparison = {
            'names': names,
            'conductivity': [],
            'surface_tension': [],
            'viscosity': [],
            'density': [],
            'M_avg': []
        }
        
        for name in names:
            liq = self.get(name)
            comparison['conductivity'].append(liq.conductivity)
            comparison['surface_tension'].append(liq.surface_tension)
            comparison['viscosity'].append(liq.viscosity)
            comparison['density'].append(liq.density)
            comparison['M_avg'].append(liq.average_molecular_weight())
        
        return comparison
    
    def find_best_for_property(self, 
                                property_name: str, 
                                maximize: bool = True) -> str:
        """Find liquid with best value for a given property"""
        values = {}
        for name, liquid in self.liquids.items():
            if hasattr(liquid, property_name):
                values[name] = getattr(liquid, property_name)
        
        if not values:
            raise ValueError(f"Property '{property_name}' not found")
        
        if maximize:
            best_name = max(values, key=values.get)
        else:
            best_name = min(values, key=values.get)
        
        return best_name


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Ionic Liquid Property Database")
    print("=" * 70)
    
    # Create database
    db = IonicLiquidDatabase()
    
    print(f"\nAvailable liquids: {', '.join(db.list_all())}")
    
    # Get specific liquid
    print("\n" + "-" * 70)
    print("EMI-Im Properties:")
    print("-" * 70)
    emi_im = db.get('EMI-Im')
    print(emi_im)
    print(f"\nAverage MW: {emi_im.average_molecular_weight()*1e3:.1f} g/mol")
    print(f"Electrochemical window: {emi_im.electrochemical_window} V")
    print(f"Thermal stability: {emi_im.thermal_stability} K")
    
    # Temperature-dependent properties
    print("\n" + "-" * 70)
    print("Temperature Dependence (EMI-Im):")
    print("-" * 70)
    for T in [273, 298, 323, 348]:
        props = emi_im.get_properties_at_T(T)
        print(f"T = {T} K:")
        print(f"  σ = {props['conductivity']:.3f} S/m")
        print(f"  μ = {props['viscosity']*1e3:.2f} mPa·s")
        print(f"  γ = {props['surface_tension']*1e3:.2f} mN/m")
        print(f"  ρ = {props['density']:.1f} kg/m³")
    
    # Comparison
    print("\n" + "-" * 70)
    print("Property Comparison:")
    print("-" * 70)
    comp = db.compare_liquids()
    print(f"{'Liquid':<12} {'σ (S/m)':<10} {'γ (mN/m)':<12} "
          f"{'μ (mPa·s)':<12} {'ρ (kg/m³)':<10}")
    print("-" * 70)
    for i, name in enumerate(comp['names']):
        print(f"{name:<12} {comp['conductivity'][i]:<10.2f} "
              f"{comp['surface_tension'][i]*1e3:<12.2f} "
              f"{comp['viscosity'][i]*1e3:<12.2f} "
              f"{comp['density'][i]:<10.0f}")
    
    # Find best liquids
    print("\n" + "-" * 70)
    print("Best Liquids for Each Property:")
    print("-" * 70)
    print(f"Highest conductivity: {db.find_best_for_property('conductivity', True)}")
    print(f"Lowest viscosity: {db.find_best_for_property('viscosity', False)}")
    print(f"Highest surface tension: {db.find_best_for_property('surface_tension', True)}")
    
    print("\n" + "=" * 70)
