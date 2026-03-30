"""
Visualization module for electrospray simulation results
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

sns.set_style("whitegrid")
sns.set_palette("husl")


def plot_voltage_sweep(results, save_path=None):
    """Plot results from voltage sweep"""
    
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    V = results['voltage']
    
    # Current vs Voltage
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(V, results['current'] * 1e9, 'o-', linewidth=2, markersize=4)
    ax1.set_xlabel('Voltage (V)', fontsize=12)
    ax1.set_ylabel('Current (nA)', fontsize=12)
    ax1.set_title('Emission Current vs Voltage', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Thrust vs Voltage
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(V, results['thrust'] * 1e6, 'o-', linewidth=2, markersize=4, color='red')
    ax2.set_xlabel('Voltage (V)', fontsize=12)
    ax2.set_ylabel('Thrust (μN)', fontsize=12)
    ax2.set_title('Thrust vs Voltage', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Specific Impulse vs Voltage
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(V, results['isp'], 'o-', linewidth=2, markersize=4, color='green')
    ax3.set_xlabel('Voltage (V)', fontsize=12)
    ax3.set_ylabel('Specific Impulse (s)', fontsize=12)
    ax3.set_title('Specific Impulse vs Voltage', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Efficiency vs Voltage
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(V, results['efficiency'] * 100, 'o-', linewidth=2, markersize=4, color='purple')
    ax4.set_xlabel('Voltage (V)', fontsize=12)
    ax4.set_ylabel('Efficiency (%)', fontsize=12)
    ax4.set_title('Efficiency vs Voltage', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim([0, 100])
    
    # Ion Fraction vs Voltage
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.plot(V, results['mode'] * 100, 'o-', linewidth=2, markersize=4, color='orange')
    ax5.set_xlabel('Voltage (V)', fontsize=12)
    ax5.set_ylabel('Ion Fraction (%)', fontsize=12)
    ax5.set_title('Ion Emission Fraction vs Voltage', fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim([0, 100])
    
    # Power vs Voltage
    ax6 = fig.add_subplot(gs[2, 1])
    power = results['current'] * results['voltage']
    ax6.plot(V, power * 1e6, 'o-', linewidth=2, markersize=4, color='brown')
    ax6.set_xlabel('Voltage (V)', fontsize=12)
    ax6.set_ylabel('Power (μW)', fontsize=12)
    ax6.set_title('Input Power vs Voltage', fontsize=13, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    
    plt.suptitle('Electrospray Performance: Voltage Sweep', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.tight_layout()
    return fig


def plot_flow_rate_sweep(results, save_path=None):
    """Plot results from flow rate sweep"""
    
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    m_dot = results['m_dot']
    
    # Current vs Flow Rate
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.loglog(m_dot, results['current'] * 1e9, 'o-', linewidth=2, markersize=4)
    ax1.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax1.set_ylabel('Current (nA)', fontsize=12)
    ax1.set_title('Emission Current vs Flow Rate', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, which='both')
    
    # Thrust vs Flow Rate
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.loglog(m_dot, results['thrust'] * 1e6, 'o-', linewidth=2, markersize=4, color='red')
    ax2.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax2.set_ylabel('Thrust (μN)', fontsize=12)
    ax2.set_title('Thrust vs Flow Rate', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')
    
    # Specific Impulse vs Flow Rate
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.semilogx(m_dot, results['isp'], 'o-', linewidth=2, markersize=4, color='green')
    ax3.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax3.set_ylabel('Specific Impulse (s)', fontsize=12)
    ax3.set_title('Specific Impulse vs Flow Rate', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Efficiency vs Flow Rate
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.semilogx(m_dot, results['efficiency'] * 100, 'o-', linewidth=2, markersize=4, color='purple')
    ax4.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax4.set_ylabel('Efficiency (%)', fontsize=12)
    ax4.set_title('Efficiency vs Flow Rate', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim([0, 100])
    
    # Jet Radius vs Flow Rate
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.loglog(m_dot, results['jet_radius'] * 1e9, 'o-', linewidth=2, markersize=4, color='orange')
    ax5.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax5.set_ylabel('Jet Radius (nm)', fontsize=12)
    ax5.set_title('Jet Radius vs Flow Rate', fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3, which='both')
    
    # Emission Mode vs Flow Rate
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.semilogx(m_dot, results['mode'] * 100, 'o-', linewidth=2, markersize=4, color='brown')
    ax6.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax6.set_ylabel('Ion Fraction (%)', fontsize=12)
    ax6.set_title('Emission Mode vs Flow Rate', fontsize=13, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim([0, 100])
    
    # Mark transition regions
    ax6.axvline(4e-12, color='red', linestyle='--', alpha=0.5, label='Pure Ion Threshold')
    ax6.axvline(1e-10, color='orange', linestyle='--', alpha=0.5, label='Droplet Transition')
    ax6.legend(fontsize=9)
    
    plt.suptitle('Electrospray Performance: Flow Rate Sweep', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.tight_layout()
    return fig


def plot_liquid_comparison(liquids_list, V_emitter=2000, m_dot=5e-11, save_path=None):
    """Compare performance of different ionic liquids"""
    
    from electrospray_simulator import SingleCapillaryEmitter
    
    results = {
        'names': [],
        'current': [],
        'thrust': [],
        'isp': [],
        'efficiency': [],
        'jet_radius': []
    }
    
    for liquid in liquids_list:
        emitter = SingleCapillaryEmitter(
            liquid=liquid,
            d_tip=20e-6,
            V_emitter=V_emitter,
            m_dot=m_dot
        )
        emitter.calculate_emission()
        
        results['names'].append(liquid.name)
        results['current'].append(emitter.I * 1e9)
        results['thrust'].append(emitter.thrust * 1e6)
        results['isp'].append(emitter.Isp)
        results['efficiency'].append(emitter.efficiency * 100)
        results['jet_radius'].append(emitter.R_jet * 1e9)
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    x_pos = np.arange(len(results['names']))
    
    # Current
    axes[0].bar(x_pos, results['current'], color='steelblue', alpha=0.8)
    axes[0].set_ylabel('Current (nA)', fontsize=11)
    axes[0].set_title('Emission Current', fontsize=12, fontweight='bold')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(results['names'], rotation=45, ha='right')
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Thrust
    axes[1].bar(x_pos, results['thrust'], color='coral', alpha=0.8)
    axes[1].set_ylabel('Thrust (μN)', fontsize=11)
    axes[1].set_title('Thrust', fontsize=12, fontweight='bold')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(results['names'], rotation=45, ha='right')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    # Specific Impulse
    axes[2].bar(x_pos, results['isp'], color='seagreen', alpha=0.8)
    axes[2].set_ylabel('Isp (s)', fontsize=11)
    axes[2].set_title('Specific Impulse', fontsize=12, fontweight='bold')
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(results['names'], rotation=45, ha='right')
    axes[2].grid(True, alpha=0.3, axis='y')
    
    # Efficiency
    axes[3].bar(x_pos, results['efficiency'], color='purple', alpha=0.8)
    axes[3].set_ylabel('Efficiency (%)', fontsize=11)
    axes[3].set_title('Propulsive Efficiency', fontsize=12, fontweight='bold')
    axes[3].set_xticks(x_pos)
    axes[3].set_xticklabels(results['names'], rotation=45, ha='right')
    axes[3].grid(True, alpha=0.3, axis='y')
    axes[3].set_ylim([0, 100])
    
    # Jet Radius
    axes[4].bar(x_pos, results['jet_radius'], color='orange', alpha=0.8)
    axes[4].set_ylabel('Jet Radius (nm)', fontsize=11)
    axes[4].set_title('Jet Radius', fontsize=12, fontweight='bold')
    axes[4].set_xticks(x_pos)
    axes[4].set_xticklabels(results['names'], rotation=45, ha='right')
    axes[4].grid(True, alpha=0.3, axis='y')
    
    # Properties comparison
    K_values = [liq.conductivity for liq in liquids_list]
    gamma_values = [liq.surface_tension * 1e3 for liq in liquids_list]
    
    ax5 = axes[5]
    x_props = np.arange(len(results['names']))
    width = 0.35
    
    ax5_twin = ax5.twinx()
    
    bars1 = ax5.bar(x_props - width/2, K_values, width, label='Conductivity', 
                    color='royalblue', alpha=0.8)
    bars2 = ax5_twin.bar(x_props + width/2, gamma_values, width, 
                         label='Surface Tension', color='crimson', alpha=0.8)
    
    ax5.set_ylabel('Conductivity (S/m)', fontsize=11, color='royalblue')
    ax5_twin.set_ylabel('Surface Tension (mN/m)', fontsize=11, color='crimson')
    ax5.set_title('Liquid Properties', fontsize=12, fontweight='bold')
    ax5.set_xticks(x_props)
    ax5.set_xticklabels(results['names'], rotation=45, ha='right')
    ax5.tick_params(axis='y', labelcolor='royalblue')
    ax5_twin.tick_params(axis='y', labelcolor='crimson')
    ax5.grid(True, alpha=0.3, axis='y')
    
    # Add legends
    lines1, labels1 = ax5.get_legend_handles_labels()
    lines2, labels2 = ax5_twin.get_legend_handles_labels()
    ax5.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)
    
    plt.suptitle(f'Ionic Liquid Performance Comparison\n' + 
                 f'(V = {V_emitter} V, ṁ = {m_dot:.1e} kg/s)',
                 fontsize=14, fontweight='bold', y=0.98)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.tight_layout()
    return fig, results


def plot_scaling_laws(liquid, save_path=None):
    """Plot theoretical scaling laws validation"""
    
    from electrospray_simulator import SingleCapillaryEmitter, ElectrosprayPhysics
    
    physics = ElectrosprayPhysics(liquid)
    
    # Flow rate range
    m_dot_range = np.logspace(-12, -9, 100)
    Q_range = m_dot_range / liquid.density
    
    # Calculate theoretical predictions
    R_jet_theory = np.array([physics.jet_radius(Q) for Q in Q_range])
    I_theory = np.array([physics.emitted_current_isothermal(Q) for Q in Q_range])
    
    # Calculate from simulation
    R_jet_sim = []
    I_sim = []
    
    for m_dot in m_dot_range:
        emitter = SingleCapillaryEmitter(
            liquid=liquid,
            m_dot=m_dot,
            V_emitter=2000
        )
        emitter.calculate_emission(use_heating=False)  # Isothermal for comparison
        R_jet_sim.append(emitter.R_jet)
        I_sim.append(emitter.I)
    
    R_jet_sim = np.array(R_jet_sim)
    I_sim = np.array(I_sim)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Jet radius scaling
    ax1 = axes[0, 0]
    ax1.loglog(m_dot_range, R_jet_theory * 1e9, 'k--', linewidth=2, label='Theory: R ∝ Q^{1/4}')
    ax1.loglog(m_dot_range, R_jet_sim * 1e9, 'ro', markersize=5, alpha=0.6, label='Simulation')
    ax1.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax1.set_ylabel('Jet Radius (nm)', fontsize=12)
    ax1.set_title('Jet Radius Scaling Law', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, which='both')
    
    # Current scaling
    ax2 = axes[0, 1]
    ax2.loglog(m_dot_range, I_theory * 1e9, 'k--', linewidth=2, label='Theory: I ∝ √ṁ')
    ax2.loglog(m_dot_range, I_sim * 1e9, 'bo', markersize=5, alpha=0.6, label='Simulation')
    ax2.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax2.set_ylabel('Current (nA)', fontsize=12)
    ax2.set_title('Current Scaling Law', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')
    
    # Droplet radius (should be ~1.89 × jet radius)
    ax3 = axes[1, 0]
    R_droplet_sim = 1.89 * R_jet_sim
    ax3.loglog(m_dot_range, R_droplet_sim * 1e6, 'go', markersize=5, alpha=0.6, label='Simulation')
    ax3.loglog(m_dot_range, 1.89 * R_jet_theory * 1e6, 'k--', linewidth=2, 
               label='Theory: R_d = 1.89R_jet')
    ax3.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax3.set_ylabel('Droplet Radius (μm)', fontsize=12)
    ax3.set_title('Droplet Size Scaling', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3, which='both')
    
    # Dimensionless parameters
    ax4 = axes[1, 1]
    Pi = np.array([physics.dimensionless_flow_rate(Q) for Q in Q_range])
    We = np.array([liquid.density * Q**2 / (2 * np.pi * R**3 * liquid.surface_tension) 
                   for Q, R in zip(Q_range, R_jet_theory)])
    
    ax4.loglog(m_dot_range, Pi, 'mo', markersize=5, alpha=0.6, label='Π (flow rate)')
    ax4.loglog(m_dot_range, We, 'co', markersize=5, alpha=0.6, label='We (Weber number)')
    ax4.axhline(1, color='red', linestyle='--', linewidth=2, alpha=0.5, label='We = 1 (stability limit)')
    ax4.set_xlabel('Mass Flow Rate (kg/s)', fontsize=12)
    ax4.set_ylabel('Dimensionless Number', fontsize=12)
    ax4.set_title('Dimensionless Parameters', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, which='both')
    
    plt.suptitle(f'Scaling Laws Validation: {liquid.name}', 
                 fontsize=15, fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.tight_layout()
    return fig


def plot_array_scaling(liquid, n_emitters_range, save_path=None):
    """Plot how array performance scales with number of emitters"""
    
    from electrospray_simulator import MultiCapillaryArray
    
    results = {
        'n_emitters': n_emitters_range,
        'thrust': [],
        'current': [],
        'power': [],
        'thrust_density': []
    }
    
    m_dot_total = 1e-8  # Fixed total flow rate
    V = 2000  # Fixed voltage
    
    for n in n_emitters_range:
        array = MultiCapillaryArray(
            liquid=liquid,
            n_emitters=n,
            V_emitter=V,
            m_dot_total=m_dot_total
        )
        array.calculate_array_emission()
        
        results['thrust'].append(array.total_thrust * 1e3)  # mN
        results['current'].append(array.total_current * 1e6)  # μA
        results['power'].append(array.total_current * V)  # W
        
        # Estimate array area (square packing)
        spacing = 1e-3  # 1 mm
        area = (spacing * np.sqrt(n))**2
        results['thrust_density'].append(array.total_thrust / area)
    
    # Convert to arrays
    for key in ['thrust', 'current', 'power', 'thrust_density']:
        results[key] = np.array(results[key])
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Thrust scaling
    axes[0, 0].plot(n_emitters_range, results['thrust'], 'o-', linewidth=2, markersize=6)
    axes[0, 0].set_xlabel('Number of Emitters', fontsize=12)
    axes[0, 0].set_ylabel('Total Thrust (mN)', fontsize=12)
    axes[0, 0].set_title('Thrust Scaling', fontsize=13, fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Current scaling
    axes[0, 1].plot(n_emitters_range, results['current'], 'o-', linewidth=2, 
                    markersize=6, color='red')
    axes[0, 1].set_xlabel('Number of Emitters', fontsize=12)
    axes[0, 1].set_ylabel('Total Current (μA)', fontsize=12)
    axes[0, 1].set_title('Current Scaling', fontsize=13, fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Power scaling
    axes[1, 0].plot(n_emitters_range, results['power'], 'o-', linewidth=2, 
                    markersize=6, color='green')
    axes[1, 0].set_xlabel('Number of Emitters', fontsize=12)
    axes[1, 0].set_ylabel('Total Power (W)', fontsize=12)
    axes[1, 0].set_title('Power Scaling', fontsize=13, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Thrust density
    axes[1, 1].plot(n_emitters_range, results['thrust_density'], 'o-', linewidth=2, 
                    markersize=6, color='purple')
    axes[1, 1].set_xlabel('Number of Emitters', fontsize=12)
    axes[1, 1].set_ylabel('Thrust Density (N/m²)', fontsize=12)
    axes[1, 1].set_title('Thrust Density', fontsize=13, fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle(f'Array Scaling: {liquid.name}\n(Total ṁ = {m_dot_total:.1e} kg/s, V = {V} V)',
                 fontsize=14, fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.tight_layout()
    return fig, results


if __name__ == "__main__":
    print("Visualization module loaded.")
    print("Available functions:")
    print("  - plot_voltage_sweep(results)")
    print("  - plot_flow_rate_sweep(results)")
    print("  - plot_liquid_comparison(liquids_list)")
    print("  - plot_scaling_laws(liquid)")
    print("  - plot_array_scaling(liquid, n_emitters_range)")
