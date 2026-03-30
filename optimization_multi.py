#!/usr/bin/env python3
"""
Multi-Parameter Optimization for Electrospray Thruster Arrays
==============================================================

Optimizes four key parameters:
1. Extractor plate spacing (gap between plates)
2. Emitter-to-extractor distance
3. Array mounting diameter (packing density)
4. Liquid propellant type

Goal: Maximize efficiency and thrust simultaneously

Author: Research analysis
Date: 2026-02-02
"""

import sys
sys.path.insert(0, '/mnt/user-data/outputs')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import itertools
from electrospray_simulator import (
    SingleCapillaryEmitter, MultiCapillaryArray,
    EMI_IM, EMI_BF4, EAN, BMI_TCM
)

plt.rcParams.update({
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'figure.dpi': 150,
})

# ============================================================
# PARAMETER SPACE DEFINITION
# ============================================================

# 1. Extractor plate spacing (assuming symmetric around emitters)
plate_spacings = np.linspace(0.3e-3, 2.0e-3, 8)  # 0.3 to 2 mm

# 2. Emitter-to-extractor distance
emitter_gaps = np.linspace(0.2e-3, 1.5e-3, 10)  # 0.2 to 1.5 mm

# 3. Array diameter (affects packing density)
array_diameters = np.linspace(5e-3, 50e-3, 8)  # 5 to 50 mm

# 4. Propellants
propellants = [EMI_IM, EMI_BF4, EAN, BMI_TCM]
propellant_names = ['EMI-Im', 'EMI-BF4', 'EAN', 'BMI-TCM']

# Fixed parameters
VOLTAGE = 2000  # V
TIP_DIAMETER = 20e-6  # 20 μm
TOTAL_FLOW = 1e-6  # 1 μg/s total
NUM_EMITTERS = 100  # Number of emitters in array

print("="*60)
print("MULTI-PARAMETER OPTIMIZATION STUDY")
print("="*60)
print(f"\nFixed Parameters:")
print(f"  Voltage: {VOLTAGE} V")
print(f"  Tip diameter: {TIP_DIAMETER*1e6:.0f} μm")
print(f"  Total flow rate: {TOTAL_FLOW*1e9:.0f} ng/s")
print(f"  Number of emitters: {NUM_EMITTERS}")
print(f"\nOptimization Variables:")
print(f"  1. Plate spacing: {plate_spacings[0]*1e3:.1f} - {plate_spacings[-1]*1e3:.1f} mm")
print(f"  2. Emitter gap: {emitter_gaps[0]*1e3:.1f} - {emitter_gaps[-1]*1e3:.1f} mm")
print(f"  3. Array diameter: {array_diameters[0]*1e3:.0f} - {array_diameters[-1]*1e3:.0f} mm")
print(f"  4. Propellants: {', '.join(propellant_names)}")

# ============================================================
# OPTIMIZATION SWEEP
# ============================================================

def calculate_packing_efficiency(array_diam, n_emitters, emitter_spacing=100e-6):
    """Calculate how efficiently emitters pack on the array"""
    array_area = np.pi * (array_diam / 2)**2
    emitter_area = n_emitters * np.pi * (emitter_spacing / 2)**2
    return min(emitter_area / array_area, 1.0)

def optimize_configuration(propellant, prop_name):
    """
    Sweep through geometric parameters for a given propellant
    Returns optimal configuration
    """
    print(f"\n{'='*60}")
    print(f"Optimizing for {prop_name}")
    print(f"{'='*60}")
    
    best_config = {
        'efficiency': 0,
        'thrust': 0,
        'isp': 0,
        'figure_of_merit': 0,
        'gap': 0,
        'plate_spacing': 0,
        'array_diam': 0,
    }
    
    results = {
        'gap': [],
        'plate_spacing': [],
        'array_diam': [],
        'efficiency': [],
        'thrust': [],
        'isp': [],
        'figure_of_merit': [],
        'packing': [],
    }
    
    total_configs = len(emitter_gaps) * len(plate_spacings) * len(array_diameters)
    count = 0
    
    for gap in emitter_gaps:
        for plate_space in plate_spacings:
            for array_diam in array_diameters:
                count += 1
                
                # Create array configuration
                try:
                    array = MultiCapillaryArray(
                        liquid=propellant,
                        n_emitters=NUM_EMITTERS,
                        d_tip=TIP_DIAMETER,
                        gap=gap,
                        V_emitter=VOLTAGE,
                        m_dot_total=TOTAL_FLOW
                    )
                    
                    # Calculate performance
                    array.calculate_array_emission(uniformity=0.95)
                    
                    efficiency = array.array_efficiency
                    thrust = array.total_thrust
                    isp = array.array_Isp
                    packing = calculate_packing_efficiency(array_diam, NUM_EMITTERS)
                    
                    # Figure of merit: weighted combination
                    # Priority: efficiency (40%), thrust (40%), compactness (20%)
                    FOM = (0.4 * efficiency + 
                           0.4 * (thrust / 1e-3) +  # Normalize to ~mN scale
                           0.2 * packing)
                    
                    # Store results
                    results['gap'].append(gap * 1e3)  # mm
                    results['plate_spacing'].append(plate_space * 1e3)
                    results['array_diam'].append(array_diam * 1e3)
                    results['efficiency'].append(efficiency * 100)
                    results['thrust'].append(thrust * 1e6)  # μN
                    results['isp'].append(isp)
                    results['figure_of_merit'].append(FOM)
                    results['packing'].append(packing * 100)
                    
                    # Update best
                    if FOM > best_config['figure_of_merit']:
                        best_config = {
                            'efficiency': efficiency * 100,
                            'thrust': thrust * 1e6,
                            'isp': isp,
                            'figure_of_merit': FOM,
                            'gap': gap * 1e3,
                            'plate_spacing': plate_space * 1e3,
                            'array_diam': array_diam * 1e3,
                            'packing': packing * 100,
                        }
                        
                except Exception as e:
                    # Skip invalid configurations
                    continue
                
                if count % 100 == 0:
                    print(f"  Progress: {count}/{total_configs} ({100*count/total_configs:.1f}%)")
    
    print(f"\n  Evaluated {len(results['gap'])} valid configurations")
    print(f"\n  OPTIMAL CONFIGURATION:")
    print(f"    Emitter-extractor gap: {best_config['gap']:.2f} mm")
    print(f"    Plate spacing: {best_config['plate_spacing']:.2f} mm")
    print(f"    Array diameter: {best_config['array_diam']:.1f} mm")
    print(f"    Packing efficiency: {best_config['packing']:.1f}%")
    print(f"\n  PERFORMANCE:")
    print(f"    Efficiency: {best_config['efficiency']:.1f}%")
    print(f"    Thrust: {best_config['thrust']:.1f} μN")
    print(f"    Isp: {best_config['isp']:.0f} s")
    print(f"    Figure of Merit: {best_config['figure_of_merit']:.3f}")
    
    return results, best_config

# Run optimization for each propellant
all_results = {}
all_best = {}

for prop, name in zip(propellants, propellant_names):
    res, best = optimize_configuration(prop, name)
    all_results[name] = res
    all_best[name] = best

# ============================================================
# VISUALIZATION
# ============================================================

print(f"\n{'='*60}")
print("Generating optimization visualizations...")
print(f"{'='*60}")

# Figure 1: Comparison of optimal configurations across propellants
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle('Optimal Configuration Comparison Across Propellants', 
             fontsize=13, fontweight='bold')

colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6']

# Geometric parameters
ax = axes[0, 0]
x = np.arange(len(propellant_names))
gaps = [all_best[name]['gap'] for name in propellant_names]
ax.bar(x, gaps, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Optimal Gap (mm)')
ax.set_title('Emitter-Extractor Gap')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.grid(axis='y', alpha=0.3)

ax = axes[0, 1]
plate_spaces = [all_best[name]['plate_spacing'] for name in propellant_names]
ax.bar(x, plate_spaces, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Optimal Spacing (mm)')
ax.set_title('Plate Spacing')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.grid(axis='y', alpha=0.3)

ax = axes[0, 2]
array_diams = [all_best[name]['array_diam'] for name in propellant_names]
ax.bar(x, array_diams, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Optimal Diameter (mm)')
ax.set_title('Array Diameter')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.grid(axis='y', alpha=0.3)

# Performance metrics
ax = axes[1, 0]
efficiencies = [all_best[name]['efficiency'] for name in propellant_names]
ax.bar(x, efficiencies, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Efficiency (%)')
ax.set_title('Thruster Efficiency')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.set_ylim(0, 100)
ax.grid(axis='y', alpha=0.3)

ax = axes[1, 1]
thrusts = [all_best[name]['thrust'] for name in propellant_names]
ax.bar(x, thrusts, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Thrust (μN)')
ax.set_title('Total Thrust')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.grid(axis='y', alpha=0.3)

ax = axes[1, 2]
isps = [all_best[name]['isp'] for name in propellant_names]
ax.bar(x, isps, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Isp (s)')
ax.set_title('Specific Impulse')
ax.set_xticks(x)
ax.set_xticklabels(propellant_names, rotation=15)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('optimization_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ Saved: optimization_comparison.png")

# Figure 2: 3D sensitivity analysis for best propellant
best_prop_name = max(all_best.keys(), key=lambda k: all_best[k]['figure_of_merit'])
print(f"\n  Best overall propellant: {best_prop_name}")

fig = plt.figure(figsize=(16, 10))
fig.suptitle(f'Sensitivity Analysis: {best_prop_name} (Best Propellant)', 
             fontsize=13, fontweight='bold')

res = all_results[best_prop_name]

# Convert to arrays for 3D plotting
gap_arr = np.array(res['gap'])
eff_arr = np.array(res['efficiency'])
thrust_arr = np.array(res['thrust'])
isp_arr = np.array(res['isp'])
fom_arr = np.array(res['figure_of_merit'])
diam_arr = np.array(res['array_diam'])

# 3D scatter: Gap vs Efficiency vs Thrust
ax = fig.add_subplot(2, 3, 1, projection='3d')
sc = ax.scatter(gap_arr, eff_arr, thrust_arr, c=fom_arr, cmap='viridis', s=20, alpha=0.6)
ax.set_xlabel('Gap (mm)')
ax.set_ylabel('Efficiency (%)')
ax.set_zlabel('Thrust (μN)')
ax.set_title('Gap vs Efficiency vs Thrust')
plt.colorbar(sc, ax=ax, label='FOM', shrink=0.5)

# 3D scatter: Gap vs Isp vs Thrust
ax = fig.add_subplot(2, 3, 2, projection='3d')
sc = ax.scatter(gap_arr, isp_arr, thrust_arr, c=fom_arr, cmap='plasma', s=20, alpha=0.6)
ax.set_xlabel('Gap (mm)')
ax.set_ylabel('Isp (s)')
ax.set_zlabel('Thrust (μN)')
ax.set_title('Gap vs Isp vs Thrust')
plt.colorbar(sc, ax=ax, label='FOM', shrink=0.5)

# 3D scatter: Array diameter vs Efficiency vs Thrust
ax = fig.add_subplot(2, 3, 3, projection='3d')
sc = ax.scatter(diam_arr, eff_arr, thrust_arr, c=fom_arr, cmap='coolwarm', s=20, alpha=0.6)
ax.set_xlabel('Array Diam (mm)')
ax.set_ylabel('Efficiency (%)')
ax.set_zlabel('Thrust (μN)')
ax.set_title('Array Size vs Performance')
plt.colorbar(sc, ax=ax, label='FOM', shrink=0.5)

# 2D contours
ax = fig.add_subplot(2, 3, 4)
# Create 2D histogram for efficiency
from scipy.stats import binned_statistic_2d
stat_eff = binned_statistic_2d(gap_arr, diam_arr, eff_arr, statistic='mean', bins=15)
im = ax.imshow(stat_eff.statistic.T, origin='lower', aspect='auto', cmap='viridis',
               extent=[gap_arr.min(), gap_arr.max(), diam_arr.min(), diam_arr.max()])
ax.set_xlabel('Gap (mm)')
ax.set_ylabel('Array Diameter (mm)')
ax.set_title('Efficiency Heatmap')
plt.colorbar(im, ax=ax, label='Efficiency (%)')

ax = fig.add_subplot(2, 3, 5)
stat_thrust = binned_statistic_2d(gap_arr, diam_arr, thrust_arr, statistic='mean', bins=15)
im = ax.imshow(stat_thrust.statistic.T, origin='lower', aspect='auto', cmap='plasma',
               extent=[gap_arr.min(), gap_arr.max(), diam_arr.min(), diam_arr.max()])
ax.set_xlabel('Gap (mm)')
ax.set_ylabel('Array Diameter (mm)')
ax.set_title('Thrust Heatmap')
plt.colorbar(im, ax=ax, label='Thrust (μN)')

ax = fig.add_subplot(2, 3, 6)
stat_fom = binned_statistic_2d(gap_arr, diam_arr, fom_arr, statistic='mean', bins=15)
im = ax.imshow(stat_fom.statistic.T, origin='lower', aspect='auto', cmap='coolwarm',
               extent=[gap_arr.min(), gap_arr.max(), diam_arr.min(), diam_arr.max()])
ax.set_xlabel('Gap (mm)')
ax.set_ylabel('Array Diameter (mm)')
ax.set_title('Figure of Merit Heatmap')
plt.colorbar(im, ax=ax, label='FOM')

# Mark optimal point
opt_gap = all_best[best_prop_name]['gap']
opt_diam = all_best[best_prop_name]['array_diam']
ax.plot(opt_gap, opt_diam, 'r*', markersize=20, markeredgecolor='white', markeredgewidth=2)
ax.text(opt_gap, opt_diam, '  OPTIMAL', color='white', fontweight='bold', fontsize=9)

plt.tight_layout()
plt.savefig('sensitivity_analysis_3d.png', dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ Saved: sensitivity_analysis_3d.png")

# ============================================================
# SUMMARY TABLE
# ============================================================

print(f"\n{'='*60}")
print("OPTIMIZATION SUMMARY")
print(f"{'='*60}")
print(f"\n{'Propellant':<12} {'Gap(mm)':<10} {'Plate(mm)':<10} {'Diam(mm)':<10} {'Eff(%)':<8} {'Thrust(μN)':<12} {'Isp(s)':<8} {'FOM':<8}")
print("-"*90)
for name in propellant_names:
    b = all_best[name]
    print(f"{name:<12} {b['gap']:<10.2f} {b['plate_spacing']:<10.2f} {b['array_diam']:<10.1f} "
          f"{b['efficiency']:<8.1f} {b['thrust']:<12.1f} {b['isp']:<8.0f} {b['figure_of_merit']:<8.3f}")

print(f"\n{'='*60}")
print(f"WINNER: {best_prop_name}")
print(f"  Recommended Configuration:")
print(f"    • Emitter-extractor gap: {all_best[best_prop_name]['gap']:.2f} mm")
print(f"    • Extractor plate spacing: {all_best[best_prop_name]['plate_spacing']:.2f} mm")
print(f"    • Array mounting diameter: {all_best[best_prop_name]['array_diam']:.1f} mm")
print(f"    • Number of emitters: {NUM_EMITTERS}")
print(f"    • Propellant: {best_prop_name}")
print(f"\n  Expected Performance:")
print(f"    • Efficiency: {all_best[best_prop_name]['efficiency']:.1f}%")
print(f"    • Total thrust: {all_best[best_prop_name]['thrust']:.1f} μN")
print(f"    • Specific impulse: {all_best[best_prop_name]['isp']:.0f} s")
print(f"    • Packing efficiency: {all_best[best_prop_name]['packing']:.1f}%")
print(f"{'='*60}")

print("\n✅ Multi-parameter optimization complete!")
