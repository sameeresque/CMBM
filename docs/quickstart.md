# CMBM Quick Start Guide

This guide will walk you through your first CMBM analysis, from installation to interpreting results.

## Prerequisites

Before starting, you should have:
- Python 3.7 or higher
- VoigtFit installed and working
- Observed spectra prepared in VoigtFit format
- Access to Cloudy photoionization model grids
- Basic familiarity with absorption line spectroscopy

## Installation

### Option 1: Install from PyPI (Recommended)

```bash
pip install cmbm
```

### Option 2: Install from Source

```bash
git clone https://github.com/yourusername/cmbm.git
cd cmbm
pip install -e .
```

### Verify Installation

```python
import cmbm
print(cmbm.__version__)
```

## Preparing Your Data

CMBM requires five input files. Here's how to prepare each:

### 1. VoigtFit Dataset (.hdf5)

Prepare your observed spectra using VoigtFit:

```python
import VoigtFit

# Create dataset
dataset = VoigtFit.DataSet(z_sys=2.5284)  # Your system redshift

# Add your observed spectra
dataset.add_data('spectrum.fits')

# Define spectral lines
dataset.add_line('HI_1215')
dataset.add_line('CII_1334')
dataset.add_line('SiII_1190')
dataset.add_line('MgII_2796')
dataset.add_line('CIV_1548')

# Save dataset
dataset.save('my_dataset.hdf5')
```

### 2. VoigtFit Results (.pickle)

Run an initial VoigtFit to get parameter estimates:

```python
# Define components for VoigtFit
dataset.add_component('SiII', z=2.52840, b=15.0, logN=13.5)
dataset.add_component('CIV', z=2.52845, b=20.0, logN=13.8)

# Fit with VoigtFit
dataset.fit()

# Save results
import pickle
with open('voigtfit_results.pickle', 'wb') as f:
    pickle.dump(dataset.best_fit, f)
```

### 3. Velocity Masks (.pkl)

Create velocity masks to exclude contaminated regions:

```python
import pandas as pd

# Define velocity ranges to FIT (not mask out)
masks = {
    'HI_1215': [[-200, -50], [50, 200]],    # Two velocity components
    'CII_1334': [[-150, 150]],               # Single broad component
    'SiII_1190': [[-100, 100]],
    'MgII_2796': [[-120, 120]],
    'CIV_1548': [[-200, 200]]
}

# Save masks
pd.to_pickle(masks, 'velocity_masks.pkl')
```

### 4. Cloudy Model Grids

You need precomputed Cloudy photoionization grids. These should be organized as:

```
/path/to/grids/
├── df_piethin/
│   ├── dfpiethin_0.2.pkl
│   ├── dfpiethin_0.2.pkl_coldens_HI
│   ├── dfpiethin_0.2.pkl_coldens_CII
│   └── ... (one file per ion)
└── df_piethick/
    ├── df_pie_thick_0.2.pkl
    ├── df_pie_thick_0.2.pkl_coldens_HI
    └── ... (one file per ion)
```

Please contact if you need additional grids.

### 5. Atomic Data File

Download the atomic transition data file (typically `atomdata.dat`) from the CMBM repository or use your own.

## Your First CMBM Fit

### Step 1: Create Configuration File

Create `my_system.yaml`:

```yaml
# Data paths
dataset_path: "my_dataset.hdf5"
voigtfit_path: "voigtfit_results.pickle"
mask_path: "velocity_masks.pkl"
transition_library_path: "atomdata.dat"

grid_paths:
  gridpath: "/path/to/your/grids/"
  gridpath_thick: "/path/to/your/grids/"
grid_z: "0.2"

# Define your system
phases:
  Low_ionization:
    - "SiII_0"      # SiII component #0
    - "MgII_1"      # MgII component #1
  High_ionization:
    - "CIV_0"       # CIV component #0

# Lines to fit
use_lines:
  - "HI_1215"
  - "CII_1334"
  - "SiII_1190"
  - "MgII_2796"
  - "CIV_1548"

# Output
output_dir: "my_first_fit"
```

### Step 2: Run the Fit

```python
import cmbm

# Create fitter
fitter = cmbm.CMBMFitter('my_system.yaml')

# Load data
fitter.load_data()

# Check setup
fitter.summary()

# Run fit (this will take 1-4 hours)
results = fitter.run_fit()
```

### Step 3: Examine Results

```python
# Get organized results
summary = fitter.get_results()

# Print evidence
print(f"Log evidence: {summary['log_evidence']:.2f} ± {summary['log_evidence_error']:.2f}")

# Print parameters
for param, stats in summary['parameter_summaries'].items():
    median = stats['median']
    lower = median - stats['percentiles']['16th']
    upper = stats['percentiles']['84th'] - median
    print(f"{param}: {median:.4f} +{upper:.4f} -{lower:.4f}")

# Get best-fit model
model_profiles, data_profiles = fitter.get_best_fit_model()
```

## Understanding the Output

CMBM creates an output directory with several files:

```
my_first_fit/
├── run.log                     # Detailed log of the fit
├── results/
│   ├── chains.txt              # MCMC chains
│   ├── params.json             # Parameter information
│   └── posterior.pdf           # Posterior distributions
└── plots/
    ├── run.pdf                 # Convergence diagnostics
    └── corner.pdf              # Corner plot
```

### Key Results to Check

1. **Log Evidence**: Bayesian evidence for model comparison
   - Higher evidence = better model fit
   - Compare different phase structures

2. **Parameter Posteriors**: Full uncertainty distributions
   - Check for multimodality
   - Look for parameter degeneracies

3. **Convergence**: Ensure the fit converged
   - Check `run.pdf` for convergence plots
   - Look for stable evidence estimate

## Understanding Parameters

Each cloud component has 5 parameters:

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| **Z** | Metallicity | -3.0 to 1.0 | [Z/H] (dex) |
| **nH** | Hydrogen density | -5.0 to 2.0 | log(cm⁻³) |
| **NHI** | HI column density | 12.0 to 22.0 | log(cm⁻²) |
| **b_turb** | Turbulent broadening | 0 to 50 | km/s |
| **z** | Redshift | z_sys ± 0.001 | - |

### Physical Interpretation

**Metallicity (Z):**
- Z = 0: Solar metallicity
- Z = -1: 0.1 solar (10% solar)
- Z = 1: 10× solar

**Hydrogen Density (nH):**
- nH = -5: Very diffuse (n ~ 10⁻⁵ cm⁻³)
- nH = 0: Moderate density (n ~ 1 cm⁻³)
- nH = 2: Dense cloud (n ~ 100 cm⁻³)

**HI Column Density (NHI):**
- NHI < 16: Optically thin
- NHI = 17-19: Intermediate
- NHI > 20: Damped Lyman-α system

## Common Workflows

### Workflow 1: Simple Two-Phase System

For a basic low/high ionization system:

```python
fitter = cmbm.CMBMFitter(
    dataset_path='data.hdf5',
    voigtfit_path='bestfit.pickle',
    mask_path='masks.pkl',
    phases={
        'Low_ion': ['SiII_0'],
        'High_ion': ['CIV_0']
    },
    use_lines=['SiII_1190', 'CIV_1548', 'HI_1215']
)
```

### Workflow 2: Multi-Component System

For multiple velocity components:

```python
phases = {
    'Component_A': ['SiII_0', 'CIV_1'],  # Blended low+high
    'Component_B': ['SiII_2', 'CIV_3'],  # Second component
}
```

### Workflow 3: Model Comparison

Compare different phase structures:

```python
# Model 1: Single phase
fitter1 = cmbm.CMBMFitter(phases={'Single': ['SiII_0']})
fitter1.load_data()
results1 = fitter1.run_fit(output_dir='model1')

# Model 2: Two phases
fitter2 = cmbm.CMBMFitter(phases={'Low': ['SiII_0'], 'High': ['CIV_1']})
fitter2.load_data()
results2 = fitter2.run_fit(output_dir='model2')

# Compare evidence
print(f"Model 1 evidence: {results1['logz']:.2f}")
print(f"Model 2 evidence: {results2['logz']:.2f}")
print(f"Bayes factor: {results2['logz'] - results1['logz']:.2f}")
```

Bayes factor interpretation:
- < 1: No preference
- 1-3: Weak evidence
- 3-5: Moderate evidence  
- \> 5: Strong evidence

## Troubleshooting

### Problem: Fit not converging

**Solution:**
- Increase `min_num_live_points` (try 800 or 1600)
- Decrease `dlogz` (try 0.3 or 0.1)
- Check your velocity masks are correct

### Problem: Parameter hitting bounds

**Solution:**
- Check VoigtFit initial values are reasonable
- Expand parameter bounds in the code
- May indicate wrong phase structure

### Problem: Very long runtime

**Solution:**
- Start with fewer lines (3-5 lines for testing)
- Reduce phase complexity (start with 2-3 components)
- Use MPI parallelization: `mpiexec -n 4 python fit_script.py`

### Problem: Poor fit quality

**Solution:**
- Check velocity masks include all relevant absorption
- Verify phase structure matches physical intuition
- Try different phase combinations
- Check for blended lines or contamination

## Best Practices

### 1. Start Simple

Begin with the simplest model:
- 2-3 phases maximum
- 5-10 spectral lines
- Well-separated velocity components

### 2. Use Good Initial Values

VoigtFit results are crucial:
- Ensure VoigtFit converged well
- Check initial parameters are physical
- Use multiple components if needed

### 3. Validate Results

Always check:
- Convergence diagnostics
- Parameter posteriors for multimodality
- Best-fit model visually matches data
- Physical reasonableness of parameters

### 4. Model Comparison

Compare models systematically:
- Start with simplest model
- Add complexity incrementally
- Use Bayesian evidence for comparison
- Consider physical interpretation

### 5. Documentation

Keep track of:
- Configuration files used
- Output directories
- Evidence values
- Physical interpretation notes

## Advanced Topics

### Parallel Execution

Speed up fitting with MPI:

```bash
mpiexec -n 8 python -c "
import cmbm
fitter = cmbm.CMBMFitter('config.yaml')
fitter.load_data()
fitter.run_fit()
"
```

### Custom Velocity Masks

Fine-tune masks for complex systems:

```python
masks = {
    'HI_1215': [
        [-300, -200],  # Component 1
        [-100, 50],    # Component 2
        [100, 250]     # Component 3
    ],
    'CIV_1548': [[-200, 200]]  # Broad single component
}
```

### High-Precision Fits

For publication-quality results:

```python
results = fitter.run_fit(
    min_num_live_points=2000,
    dlogz=0.1,
    min_ess=1000
)
```

## Next Steps

- Read the full [API documentation](api.md)
- Check out [advanced examples](../examples/)
- Join the [GitHub discussions](https://github.com/yourusername/cmbm/discussions)
- Report issues on [GitHub](https://github.com/yourusername/cmbm/issues)

## Getting Help

- **Documentation**: [https://cmbm.readthedocs.io/](https://cmbm.readthedocs.io/)
- **GitHub Issues**: Report bugs and request features
- **Discussions**: Ask questions and share results
- **Email**: contact@yourdomain.com (for private inquiries)

Happy fitting! 🌟
