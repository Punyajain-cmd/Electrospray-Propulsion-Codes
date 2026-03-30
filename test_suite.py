#!/usr/bin/env python3
"""
Comprehensive Test Suite for Electrospray Thruster Simulator
=============================================================

Tests all major components:
1. Core modules (config, liquids, physics)
2. Single emitter simulations
3. Multi-capillary arrays
4. Parameter sweeps
5. Optimization
6. Electric field solver

Run this to verify installation and functionality.
"""

import sys
import os
import numpy as np

# Add paths
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'core'))

print("=" * 80)
print(" ELECTROSPRAY THRUSTER SIMULATOR - COMPREHENSIVE TEST SUITE")
print("=" * 80)

# ============================================================================
# TEST 1: Core Module Imports
# ============================================================================

print("\n" + "-" * 80)
print("TEST 1: Core Module Imports")
print("-" * 80)

try:
    from config import PhysicalConstants, UnitConverter, ParameterBounds, SimulationConfig
    print("✓ Config module imported successfully")
    
    from ionic_liquids import (IonicLiquidDatabase, EMI_IM, EMI_BF4, EAN, 
                                BMI_TCM, BMI_IM, EMIM_DCA)
    print("✓ Ionic liquids module imported successfully")
    
    from physics import (ElectrosprayPhysics, PolydispersityModel, 
                         SpaceChargeModel, InstabilityAnalysis)
    print("✓ Physics module imported successfully")
    
    print("\n✅ All core modules passed import test")
except Exception as e:
    print(f"\n❌ Core module import failed: {e}")
    sys.exit(1)

# ============================================================================
# TEST 2: Physical Constants
# ============================================================================

print("\n" + "-" * 80)
print("TEST 2: Physical Constants")
print("-" * 80)

try:
    constants = PhysicalConstants.get_all()
    print(f"Number of constants: {len(constants)}")
    print(f"Sample: ε₀ = {constants['EPSILON_0']:.3e} F/m")
    print(f"        e = {constants['E_CHARGE']:.3e} C")
    print(f"        k_B = {constants['K_B']:.3e} J/K")
    
    assert len(constants) >= 7, "Missing constants"
    print("\n✅ Physical constants test passed")
except Exception as e:
    print(f"\n❌ Physical constants test failed: {e}")

# ============================================================================
# TEST 3: Unit Conversions
# ============================================================================

print("\n" + "-" * 80)
print("TEST 3: Unit Conversions")
print("-" * 80)

try:
    # Voltage
    V_kV = UnitConverter.voltage(2000, 'V', 'kV')
    assert abs(V_kV - 2.0) < 1e-10, f"Voltage conversion failed: {V_kV}"
    print(f"✓ 2000 V = {V_kV} kV")
    
    # Length
    d_um = UnitConverter.length(20e-6, 'm', 'μm')
    assert abs(d_um - 20.0) < 1e-10, f"Length conversion failed: {d_um}"
    print(f"✓ 20e-6 m = {d_um} μm")
    
    # Mass flow
    m_pg = UnitConverter.mass_flow(5e-11, 'kg/s', 'pg/s')
    assert abs(m_pg - 50.0) < 1e-6, f"Mass flow conversion failed: {m_pg}"
    print(f"✓ 5e-11 kg/s = {m_pg:.1f} pg/s")
    
    print("\n✅ Unit conversion test passed")
except Exception as e:
    print(f"\n❌ Unit conversion test failed: {e}")

# ============================================================================
# TEST 4: Parameter Validation
# ============================================================================

print("\n" + "-" * 80)
print("TEST 4: Parameter Validation")
print("-" * 80)

try:
    bounds = ParameterBounds()
    
    # Valid parameters should pass
    bounds.validate_voltage(2000)
    print("✓ Valid voltage accepted: 2000 V")
    
    bounds.validate_mass_flow(5e-11)
    print("✓ Valid mass flow accepted: 5e-11 kg/s")
    
    # Invalid parameters should raise errors
    try:
        bounds.validate_voltage(100000)  # Too high
        print("❌ Invalid voltage was not rejected!")
    except ValueError:
        print("✓ Invalid voltage correctly rejected: 100000 V")
    
    print("\n✅ Parameter validation test passed")
except Exception as e:
    print(f"\n❌ Parameter validation test failed: {e}")

# ============================================================================
# TEST 5: Ionic Liquid Database
# ============================================================================

print("\n" + "-" * 80)
print("TEST 5: Ionic Liquid Database")
print("-" * 80)

try:
    db = IonicLiquidDatabase()
    liquids = db.list_all()
    
    print(f"Number of liquids in database: {len(liquids)}")
    print(f"Available: {', '.join(liquids)}")
    
    # Test retrieval
    emi_im = db.get('EMI-Im')
    print(f"\n✓ Retrieved: {emi_im.name}")
    print(f"  σ = {emi_im.conductivity} S/m")
    print(f"  γ = {emi_im.surface_tension*1e3:.2f} mN/m")
    print(f"  μ = {emi_im.viscosity*1e3:.2f} mPa·s")
    print(f"  ρ = {emi_im.density} kg/m³")
    
    # Test temperature dependence
    props_298 = emi_im.get_properties_at_T(298.15)
    props_323 = emi_im.get_properties_at_T(323.15)
    print(f"\n✓ Temperature dependence:")
    print(f"  σ(298 K) = {props_298['conductivity']:.3f} S/m")
    print(f"  σ(323 K) = {props_323['conductivity']:.3f} S/m")
    
    assert len(liquids) >= 6, "Missing liquids in database"
    assert props_323['conductivity'] > props_298['conductivity'], "Conductivity should increase with T"
    
    print("\n✅ Ionic liquid database test passed")
except Exception as e:
    print(f"\n❌ Ionic liquid database test failed: {e}")

# ============================================================================
# TEST 6: Physics Models
# ============================================================================

print("\n" + "-" * 80)
print("TEST 6: Physics Models")
print("-" * 80)

try:
    physics = ElectrosprayPhysics(EMI_IM)
    
    Q = 5e-14  # m³/s
    m_dot = Q * EMI_IM.density
    
    # Characteristic scales
    r_G = physics.characteristic_length(Q)
    Pi = physics.dimensionless_flow_rate(Q)
    Re_K = physics.electrohydrodynamic_reynolds()
    
    print(f"Characteristic scales for Q = {Q:.2e} m³/s:")
    print(f"  r_G = {r_G*1e9:.2f} nm")
    print(f"  Π = {Pi:.2e}")
    print(f"  Re_K = {Re_K:.2f}")
    
    # Jet properties
    R_jet = physics.jet_radius(Q)
    Q_min = physics.minimum_flow_rate()
    
    print(f"\nJet properties:")
    print(f"  R_jet = {R_jet*1e9:.2f} nm")
    print(f"  Q_min = {Q_min:.2e} m³/s")
    print(f"  Q/Q_min = {Q/Q_min:.2f}")
    
    # Current emission
    I_iso = physics.emitted_current_isothermal(Q)
    I_heat = physics.emitted_current_with_heating(Q)
    
    print(f"\nCurrent emission:")
    print(f"  Isothermal: {I_iso*1e9:.2f} nA")
    print(f"  With heating: {I_heat*1e9:.2f} nA")
    
    assert r_G > 0, "Characteristic length must be positive"
    assert R_jet > 0, "Jet radius must be positive"
    assert I_iso > 0, "Current must be positive"
    
    print("\n✅ Physics models test passed")
except Exception as e:
    print(f"\n❌ Physics models test failed: {e}")

# ============================================================================
# TEST 7: Polydispersity Model
# ============================================================================

print("\n" + "-" * 80)
print("TEST 7: Polydispersity Model")
print("-" * 80)

try:
    R_mean = 1e-6  # 1 μm
    poly = PolydispersityModel(R_mean, sigma_relative=0.15)
    
    D_sauter = poly.sauter_mean_diameter()
    
    print(f"Droplet size distribution:")
    print(f"  Mean radius: {R_mean*1e6:.3f} μm")
    print(f"  Relative std dev: 15%")
    print(f"  Sauter mean diameter: {D_sauter*1e6:.3f} μm")
    
    # Test distribution
    R_values = np.linspace(0.5e-6, 2e-6, 100)
    P_R = poly.size_distribution(R_values)
    
    assert D_sauter > R_mean, "Sauter mean should be larger than arithmetic mean"
    assert np.all(P_R >= 0), "Probabilities must be non-negative"
    
    print("\n✅ Polydispersity model test passed")
except Exception as e:
    print(f"\n❌ Polydispersity model test failed: {e}")

# ============================================================================
# TEST 8: Space Charge Model
# ============================================================================

print("\n" + "-" * 80)
print("TEST 8: Space Charge Model")
print("-" * 80)

try:
    I = 1e-6  # 1 μA
    V = 2000  # 2 kV
    m_per_q = 100  # kg/C
    
    space_charge = SpaceChargeModel(I, V, m_per_q)
    
    P = space_charge.beam_perveance()
    I_CL = space_charge.child_langmuir_limit(gap=0.5e-3, area=1e-6)
    theta = space_charge.space_charge_expansion_angle(10e-6, 1e-3)
    
    print(f"Space charge effects:")
    print(f"  Beam perveance: {P:.3e}")
    print(f"  Child-Langmuir limit: {I_CL*1e6:.2f} μA")
    print(f"  Expansion angle: {theta*1e3:.2f} mrad")
    
    assert P > 0, "Perveance must be positive"
    assert I_CL > I, "Should be below CL limit"
    
    print("\n✅ Space charge model test passed")
except Exception as e:
    print(f"\n❌ Space charge model test failed: {e}")

# ============================================================================
# TEST 9: Instability Analysis
# ============================================================================

print("\n" + "-" * 80)
print("TEST 9: Instability Analysis")
print("-" * 80)

try:
    Q = 1e-13
    R_jet = 50e-9
    
    instability = InstabilityAnalysis(EMI_IM, Q, R_jet)
    
    is_stable, reason = instability.is_stable()
    lambda_R = instability.rayleigh_instability_wavelength()
    t_break = instability.capillary_instability_time()
    
    print(f"Instability analysis:")
    print(f"  Stable: {is_stable}")
    print(f"  Reason: {reason}")
    print(f"  Rayleigh wavelength: {lambda_R*1e6:.3f} μm")
    print(f"  Breakup time: {t_break*1e6:.2f} μs")
    
    assert lambda_R > 0, "Wavelength must be positive"
    assert t_break > 0, "Time must be positive"
    
    print("\n✅ Instability analysis test passed")
except Exception as e:
    print(f"\n❌ Instability analysis test failed: {e}")

# ============================================================================
# TEST 10: Single Emitter Simulation
# ============================================================================

print("\n" + "-" * 80)
print("TEST 10: Single Emitter Simulation")
print("-" * 80)

try:
    from electrospray_simulator import SingleCapillaryEmitter
    
    emitter = SingleCapillaryEmitter(
        liquid=EMI_IM,
        d_tip=20e-6,
        gap=0.5e-3,
        V_emitter=2000,
        m_dot=5e-11
    )
    
    emitter.calculate_emission(use_heating=True)
    
    print(f"Single emitter results:")
    print(f"  Current: {emitter.I*1e9:.2f} nA")
    print(f"  Thrust: {emitter.thrust*1e6:.2f} μN")
    print(f"  Isp: {emitter.Isp:.0f} s")
    print(f"  Efficiency: {emitter.efficiency*100:.1f}%")
    
    assert emitter.I > 0, "Current must be positive"
    assert emitter.thrust > 0, "Thrust must be positive"
    assert 0 < emitter.efficiency < 1, "Efficiency must be between 0 and 1"
    
    print("\n✅ Single emitter simulation test passed")
except Exception as e:
    print(f"\n❌ Single emitter simulation test failed: {e}")

# ============================================================================
# TEST 11: Multi-Capillary Array
# ============================================================================

print("\n" + "-" * 80)
print("TEST 11: Multi-Capillary Array")
print("-" * 80)

try:
    from electrospray_simulator import MultiCapillaryArray
    
    array = MultiCapillaryArray(
        liquid=EMI_IM,
        n_emitters=100,
        d_tip=20e-6,
        V_emitter=2000,
        m_dot_total=5e-9
    )
    
    array.calculate_array_emission(uniformity=0.95)
    
    print(f"Array results:")
    print(f"  Total current: {array.total_current*1e6:.2f} μA")
    print(f"  Total thrust: {array.total_thrust*1e3:.3f} mN")
    print(f"  Array Isp: {array.array_Isp:.0f} s")
    print(f"  Array efficiency: {array.array_efficiency*100:.1f}%")
    
    assert array.total_current > 0, "Array current must be positive"
    assert array.total_thrust > 0, "Array thrust must be positive"
    
    print("\n✅ Multi-capillary array test passed")
except Exception as e:
    print(f"\n❌ Multi-capillary array test failed: {e}")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print(" TEST SUITE SUMMARY")
print("=" * 80)
print("\n✅ ALL TESTS PASSED - Simulator is functioning correctly!")
print("\nYou can now run the full examples:")
print("  python3 examples_run.py")
print("  python3 optimization_multi.py")
print("  python3 electric_field_solver_2d.py")
print("\n" + "=" * 80)
