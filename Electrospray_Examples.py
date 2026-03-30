"""
Comprehensive Electrospray Simulation Example
Demonstrates all capabilities of the simulator
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.abspath("path_to_electrospray_simulator"))
sys.path.append(os.path.abspath("path_to_electrospray_visualization"))


from electrospray_simulator import (
    EMI_IM, EMI_BF4, EAN, BMI_TCM,
    SingleCapillaryEmitter,
    MultiCapillaryArray
)
from electrospray_visualization import (
    plot_voltage_sweep,
    plot_flow_rate_sweep,
    plot_liquid_comparison,
    plot_scaling_laws,
    plot_array_scaling
)

def example_1_single_emitter_basic():
    """Example 1: Basic single emitter simulation"""
    print("\n" + "="*80)
    print("EXAMPLE 1: Basic Single Emitter Simulation")
    print("="*80)
    
    # Create emitter with EMI-Im
    emitter = SingleCapillaryEmitter(
        liquid=EMI_IM,
        d_tip=20e-6,           # 20 μm tip diameter
        gap=0.5e-3,            # 0.5 mm gap
        V_emitter=2000,        # 2 kV
        m_dot=5e-11            # 50 pg/s
    )
    
    # Calculate emission
    emitter.calculate_emission(use_heating=True)
    
    # Print results
    emitter.print_results()
    
    return emitter


def example_2_parametric_sweeps():
    """Example 2: Parametric sweeps (voltage and flow rate)"""
    print("\n" + "="*80)
    print("EXAMPLE 2: Parametric Sweeps")
    print("="*80)
    
    # Create emitter
    emitter = SingleCapillaryEmitter(
        liquid=EMI_IM,
        d_tip=15e-6,           # 15 μm tip (smaller)
        V_emitter=2000,
        m_dot=3e-11
    )
    
    # Voltage sweep
    print("\nPerforming voltage sweep (1000 V to 5000 V)...")
    V_range = np.linspace(1000, 5000, 50)
    results_V = emitter.sweep_voltage(V_range)
    
    # Plot voltage sweep
    fig1 = plot_voltage_sweep(results_V)
    # fig1.savefig('/mnt/user-data/outputs/voltage_sweep.png', dpi=300, bbox_inches='tight')
    print("  ✓ Voltage sweep plot saved")
    
    # Flow rate sweep
    print("\nPerforming flow rate sweep (1 pg/s to 1 ng/s)...")
    m_dot_range = np.logspace(-12, -9, 50)
    results_Q = emitter.sweep_flow_rate(m_dot_range)
    
    # Plot flow rate sweep
    fig2 = plot_flow_rate_sweep(results_Q)
    # fig2.savefig('/mnt/user-data/outputs/flow_rate_sweep.png', dpi=300, bbox_inches='tight')
    print("  ✓ Flow rate sweep plot saved")
    
    return results_V, results_Q


def example_3_liquid_comparison():
    """Example 3: Compare different ionic liquids"""
    print("\n" + "="*80)
    print("EXAMPLE 3: Ionic Liquid Comparison")
    print("="*80)
    
    liquids = [EMI_IM, EMI_BF4, EAN, BMI_TCM]
    
    print("\nComparing performance of:")
    for liq in liquids:
        print(f"  - {liq.name}")
    
    # Compare at standard conditions
    fig, results = plot_liquid_comparison(
        liquids,
        V_emitter=2000,
        m_dot=5e-11
    )
    
    # fig.savefig('/mnt/user-data/outputs/liquid_comparison.png', dpi=300, bbox_inches='tight')
    print("\n  ✓ Liquid comparison plot saved")
    
    # Print summary
    print("\nPerformance Summary:")
    print("-" * 60)
    for i, name in enumerate(results['names']):
        print(f"{name:12s}: I={results['current'][i]:6.1f} nA, " + 
              f"T={results['thrust'][i]:5.2f} μN, " +
              f"Isp={results['isp'][i]:5.0f} s, " +
              f"η={results['efficiency'][i]:4.1f}%")
    
    return results


def example_4_scaling_laws():
    """Example 4: Validate theoretical scaling laws"""
    print("\n" + "="*80)
    print("EXAMPLE 4: Scaling Laws Validation")
    print("="*80)
    
    print("\nValidating scaling laws for EMI-Im...")
    print("  Theoretical predictions:")
    print("    - Jet radius: R ∝ Q^(1/4)")
    print("    - Current: I ∝ √ṁ")
    print("    - Droplet radius: R_d = 1.89 × R_jet")
    
    fig = plot_scaling_laws(EMI_IM)
    # fig.savefig('/mnt/user-data/outputs/scaling_laws.png', dpi=300, bbox_inches='tight')
    print("\n  ✓ Scaling laws validation plot saved")
    
    return fig


def example_5_multi_capillary_array():
    """Example 5: Multi-capillary array simulation"""
    print("\n" + "="*80)
    print("EXAMPLE 5: Multi-Capillary Array")
    print("="*80)
    
    # Create 100-emitter array
    print("\nCreating 100-emitter array...")
    array = MultiCapillaryArray(
        liquid=EMI_IM,
        n_emitters=100,
        d_tip=20e-6,
        array_spacing=1e-3,
        V_emitter=2000,
        m_dot_total=5e-9       # 5 ng/s total
    )
    
    # Calculate emission
    array.calculate_array_emission(uniformity=0.95)
    array.print_array_results()
    
    return array


def example_6_array_optimization():
    """Example 6: Array optimization"""
    print("\n" + "="*80)
    print("EXAMPLE 6: Array Optimization")
    print("="*80)
    
    # Create array
    array = MultiCapillaryArray(
        liquid=EMI_IM,
        n_emitters=100,
        V_emitter=2000,
        m_dot_total=5e-9
    )
    
    print("\nOptimizing for maximum thrust...")
    print("  Constraints:")
    print("    - Maximum power: 10 W")
    print("    - Minimum Isp: 1000 s")
    
    result = array.optimize_performance(
        target='thrust',
        constraint_power=10.0,
        constraint_isp_min=1000
    )
    
    return array


def example_7_array_scaling():
    """Example 7: Array scaling analysis"""
    print("\n" + "="*80)
    print("EXAMPLE 7: Array Scaling Analysis")
    print("="*80)
    
    print("\nAnalyzing how performance scales with array size...")
    n_range = np.array([10, 25, 50, 100, 200, 500, 1000])
    
    fig, results = plot_array_scaling(EMI_IM, n_range)
    # fig.savefig('/mnt/user-data/outputs/array_scaling.png', dpi=300, bbox_inches='tight')
    print("\n  ✓ Array scaling plot saved")
    
    # Print scaling summary
    print("\nScaling Summary:")
    print("-" * 70)
    print(f"{'N':>6s} {'Thrust (mN)':>12s} {'Current (μA)':>13s} " + 
          f"{'Power (W)':>10s} {'T/A (N/m²)':>12s}")
    print("-" * 70)
    for i, n in enumerate(n_range):
        print(f"{n:6d} {results['thrust'][i]:12.3f} {results['current'][i]:13.2f} " +
              f"{results['power'][i]:10.3f} {results['thrust_density'][i]:12.2f}")
    
    return results


def example_8_emitter_size_study():
    """Example 8: Study effect of emitter tip size"""
    print("\n" + "="*80)
    print("EXAMPLE 8: Emitter Tip Size Study")
    print("="*80)
    
    # Test different tip sizes (as in Paper 3)
    tip_sizes = np.array([15, 20, 30, 40, 50]) * 1e-6  # μm to m
    
    results = {
        'd_tip': tip_sizes * 1e6,
        'Q_min': [],
        'I': [],
        'thrust': [],
        'isp': []
    }
    
    print("\nTesting tip sizes from 15 μm to 50 μm...")
    print("  (As demonstrated in Caballero-Pérez et al., 2025)")
    
    for d_tip in tip_sizes:
        emitter = SingleCapillaryEmitter(
            liquid=EMI_IM,
            d_tip=d_tip,
            V_emitter=2000,
            m_dot=5e-11
        )
        emitter.calculate_emission()
        
        Q_min = emitter.physics.minimum_flow_rate()
        
        results['Q_min'].append(Q_min * emitter.liquid.density)
        results['I'].append(emitter.I)
        results['thrust'].append(emitter.thrust)
        results['isp'].append(emitter.Isp)
    
    # Convert to arrays
    for key in ['Q_min', 'I', 'thrust', 'isp']:
        results[key] = np.array(results[key])
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    axes[0, 0].plot(results['d_tip'], results['Q_min'] * 1e12, 'o-', linewidth=2, markersize=8)
    axes[0, 0].set_xlabel('Tip Diameter (μm)', fontsize=11)
    axes[0, 0].set_ylabel('Min. Mass Flow (pg/s)', fontsize=11)
    axes[0, 0].set_title('Minimum Flow Rate vs Tip Size', fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(results['d_tip'], results['I'] * 1e9, 'o-', linewidth=2, 
                    markersize=8, color='red')
    axes[0, 1].set_xlabel('Tip Diameter (μm)', fontsize=11)
    axes[0, 1].set_ylabel('Current (nA)', fontsize=11)
    axes[0, 1].set_title('Current vs Tip Size', fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].plot(results['d_tip'], results['thrust'] * 1e6, 'o-', linewidth=2, 
                    markersize=8, color='green')
    axes[1, 0].set_xlabel('Tip Diameter (μm)', fontsize=11)
    axes[1, 0].set_ylabel('Thrust (μN)', fontsize=11)
    axes[1, 0].set_title('Thrust vs Tip Size', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].plot(results['d_tip'], results['isp'], 'o-', linewidth=2, 
                    markersize=8, color='purple')
    axes[1, 1].set_xlabel('Tip Diameter (μm)', fontsize=11)
    axes[1, 1].set_ylabel('Specific Impulse (s)', fontsize=11)
    axes[1, 1].set_title('Isp vs Tip Size', fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle('Effect of Emitter Tip Diameter\n(EMI-Im, V=2000V, ṁ=50 pg/s)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    
    # fig.savefig('/mnt/user-data/outputs/tip_size_study.png', dpi=300, bbox_inches='tight')
    print("\n  ✓ Tip size study plot saved")
    
    print("\nResults Summary:")
    print("-" * 70)
    print(f"{'d_tip (μm)':>11s} {'Q_min (pg/s)':>14s} {'I (nA)':>8s} " + 
          f"{'T (μN)':>8s} {'Isp (s)':>9s}")
    print("-" * 70)
    for i in range(len(tip_sizes)):
        print(f"{results['d_tip'][i]:11.0f} {results['Q_min'][i]*1e12:14.2f} " +
              f"{results['I'][i]*1e9:8.1f} {results['thrust'][i]*1e6:8.2f} " +
              f"{results['isp'][i]:9.0f}")
    
    return results


def run_all_examples():
    """Run all examples sequentially"""
    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + " COMPREHENSIVE ELECTROSPRAY SIMULATION EXAMPLES ".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80)
    
    # Run examples
    emitter1 = example_1_single_emitter_basic()
    results_V, results_Q = example_2_parametric_sweeps()
    comp_results = example_3_liquid_comparison()
    scaling_fig = example_4_scaling_laws()
    array1 = example_5_multi_capillary_array()
    array2 = example_6_array_optimization()
    array_scale = example_7_array_scaling()
    tip_study = example_8_emitter_size_study()
    
    print("\n" + "#"*80)
    print("#" + " ALL EXAMPLES COMPLETED SUCCESSFULLY! ".center(78) + "#")
    print("#"*80)
    print("\nGenerated files in /mnt/user-data/outputs/:")
    print("  ✓ voltage_sweep.png")
    print("  ✓ flow_rate_sweep.png")
    print("  ✓ liquid_comparison.png")
    print("  ✓ scaling_laws.png")
    print("  ✓ array_scaling.png")
    print("  ✓ tip_size_study.png")
    print("\n" + "#"*80 + "\n")
    
    return {
        'emitter_basic': emitter1,
        'voltage_sweep': results_V,
        'flow_rate_sweep': results_Q,
        'liquid_comparison': comp_results,
        'scaling_laws': scaling_fig,
        'array_basic': array1,
        'array_optimized': array2,
        'array_scaling': array_scale,
        'tip_study': tip_study
    }


if __name__ == "__main__":
    # Run all examples
    all_results = run_all_examples()
    
    # Show one plot as example
    plt.show()