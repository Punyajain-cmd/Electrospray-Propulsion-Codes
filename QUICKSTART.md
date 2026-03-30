# Quick Start Guide
## Electrospray Thruster Simulation Suite v2.0

Get up and running in 5 minutes!

---

## 🚀 Step 1: Installation (1 minute)

```bash
# Extract the suite
cd electrospray_thruster_suite/

# Install dependencies
pip install numpy scipy matplotlib seaborn --break-system-packages

# Test installation
python3 test_suite.py
```

Expected output:
```
✅ ALL TESTS PASSED - Simulator is functioning correctly!
```

---

## 📚 Step 2: Your First Simulation (2 minutes)

Create a file `my_first_simulation.py`:

```python
# Import modules
from core.ionic_liquids import IonicLiquidDatabase
from electrospray_simulator import SingleCapillaryEmitter

# Get ionic liquid
db = IonicLiquidDatabase()
liquid = db.get('EMI-Im')  # Standard reference liquid

# Create emitter
emitter = SingleCapillaryEmitter(
    liquid=liquid,
    d_tip=20e-6,       # 20 μm tip diameter
    gap=0.5e-3,        # 0.5 mm emitter-extractor gap
    V_emitter=2000,    # 2 kV applied voltage
    m_dot=5e-11        # 50 pg/s mass flow rate
)

# Calculate emission
emitter.calculate_emission(use_heating=True)

# View results
emitter.print_results()
```

Run it:
```bash
python3 my_first_simulation.py
```

You'll see complete output including current, thrust, Isp, and efficiency!

---

## 📊 Step 3: Generate Plots (2 minutes)

Add to your script:

```python
import numpy as np
from electrospray_visualization import plot_voltage_sweep
import matplotlib.pyplot as plt

# Sweep voltage from 1 kV to 5 kV
V_range = np.linspace(1000, 5000, 50)
results = emitter.sweep_voltage(V_range)

# Create beautiful plots
fig = plot_voltage_sweep(results)
fig.savefig('my_voltage_sweep.png', dpi=300, bbox_inches='tight')
print("Plot saved: my_voltage_sweep.png")

plt.show()
```

This generates 6 plots showing how current, thrust, Isp, efficiency, and power vary with voltage!

---

## 🎯 What's Next?

### Try Different Liquids
```python
# Compare all 6 available liquids
for name in db.list_all():
    liquid = db.get(name)
    emitter = SingleCapillaryEmitter(liquid=liquid, V_emitter=2000, m_dot=5e-11)
    emitter.calculate_emission()
    print(f"{name:10s}: I={emitter.I*1e9:6.1f} nA, T={emitter.thrust*1e6:5.2f} μN")
```

### Run Multi-Capillary Array
```python
from electrospray_simulator import MultiCapillaryArray

array = MultiCapillaryArray(
    liquid=liquid,
    n_emitters=100,      # 100-emitter array
    V_emitter=2000,
    m_dot_total=5e-9     # 5 ng/s total
)

array.calculate_array_emission(uniformity=0.95)
array.print_array_results()
```

### Use Advanced Features
```python
from core.physics import PolydispersityModel, SpaceChargeModel, InstabilityAnalysis

# Droplet size distribution
poly = PolydispersityModel(R_mean=1e-6, sigma_relative=0.15)
D_sauter = poly.sauter_mean_diameter()

# Space charge effects
space_charge = SpaceChargeModel(I=1e-6, V=2000, m_per_q=100)
I_limit = space_charge.child_langmuir_limit(gap=0.5e-3, area=1e-6)

# Stability check
instability = InstabilityAnalysis(liquid, Q=1e-13, R_jet=50e-9)
is_stable, reason = instability.is_stable()
print(f"Stable: {is_stable} - {reason}")
```

### Run Complete Examples
```bash
# Run all 8 comprehensive examples
python3 examples_run.py

# Run multi-parameter optimization
python3 optimization_multi.py

# Run 2D electric field solver
python3 electric_field_solver_2d.py
```

---

## 📖 Learn More

- **Full Documentation**: See `README.md`
- **All Features**: See `IMPROVEMENTS.md`
- **Examples**: Look at `examples_run.py`
- **API Reference**: Check docstrings in each module

---

## 🆘 Common Issues

**Problem:** Import errors  
**Solution:** Make sure you're in the `electrospray_thruster_suite/` directory

**Problem:** "Module not found"  
**Solution:** Check Python path:
```python
import sys
sys.path.insert(0, 'core')
```

**Problem:** Plots don't show  
**Solution:** Add `plt.show()` at the end

**Problem:** Warnings about unstable flow  
**Solution:** This is normal for very low flow rates. Increase `m_dot` or check `Q_min`

---

## ✅ You're Ready!

You now know how to:
- ✓ Run basic simulations
- ✓ Generate plots
- ✓ Compare liquids
- ✓ Use arrays
- ✓ Access advanced features

**Happy simulating!** 🚀

---

**Need help?** Check the comprehensive `README.md` or examine the test suite examples in `test_suite.py`.
