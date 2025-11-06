# CMBM Documentation

**Cloud-by-cloud Multiphase Bayesian ionization Modeling**

Welcome to the CMBM documentation! CMBM is a Python package for Bayesian analysis of astronomical absorption line systems using nested sampling and Cloudy photoionization models.

## What is CMBM?

CMBM provides a framework for modeling absorption systems in quasar spectra, circumgalactic medium observations, and interstellar absorption studies. It combines:

- **Physical modeling**: Cloudy photoionization grids for column density predictions
- **Bayesian inference**: UltraNest nested sampling for robust parameter estimation
- **Multi-phase approach**: Simultaneous fitting of multiple ionization phases
- **Cloud-by-cloud methodology**: Individual velocity components with distinct physical conditions

### Key Features

✨ **Multi-phase Modeling**: Fit low, intermediate, and high ionization phases simultaneously

🎯 **Robust Uncertainties**: Full posterior distributions via nested sampling

🔬 **Physical Parameters**: Metallicity, density, column density, turbulence, and redshift

📊 **Model Comparison**: Bayesian evidence for objective model selection

⚡ **Parallel Execution**: MPI support for efficient computation

🔧 **Flexible Configuration**: YAML-based setup or direct Python API

## Quick Links

- **[Installation Guide](installation.md)** - Get CMBM up and running
- **[Quick Start Tutorial](quickstart.md)** - Your first CMBM analysis
- **[API Reference](api.md)** - Complete function and class documentation
- **[Examples](../examples/)** - Working code examples
- **[GitHub Repository](https://github.com/yourusername/cmbm)** - Source code and issues

## Getting Started

### Installation

```bash
pip install cmbm
```

### Your First Fit

```python
import cmbm

# Create fitter from configuration file
fitter = cmbm.CMBMFitter('my_config.yaml')

# Load data
fitter.load_data()

# Run nested sampling
results = fitter.run_fit()

# Get results
summary = fitter.get_results()
fitter.summary()
```

See the [Quick Start Guide](quickstart.md) for a complete tutorial.

## Methodology
### Physical Parameters

Each cloud component in CMBM is characterized by 5 parameters:

| Parameter | Symbol | Description | Physical Range |
|-----------|--------|-------------|----------------|
| Metallicity | Z | Metal abundance relative to solar | -4 to +1.5 [Z/H] |
| Hydrogen Density | nH | Volume density of hydrogen | -6 to +2 log(cm⁻³) |
| HI Column Density | NHI | Neutral hydrogen column | 12 to 23 log(cm⁻²) |
| Turbulent Velocity | b_turb | Non-thermal broadening | 0 to 200 km/s |
| Redshift | z | Line-of-sight velocity | z_sys ± Δz |

### Cloudy Grid Interpolation

CMBM uses precomputed Cloudy photoionization grids to predict column densities:

1. **Input**: Z, nH, NHI for a cloud
2. **Cloudy**: Solves ionization equilibrium
3. **Output**: Column densities for all ions of interest (CII, CIV, SiII, etc.)
4. **Interpolation**: Fast lookup during nested sampling

This ensures physical self-consistency between different ions in the same cloud.

## Example

### 1. Circumgalactic Medium (CGM) Studies

Model the multi-phase structure of CGM absorbers around galaxies:

```python
fitter = cmbm.CMBMFitter(
    phases={
        'Cool_phase': ['MgII_0', 'SiII_1', 'FeII_2'],
        'Warm_phase': ['CII_0', 'SiIII_1'],
        'Hot_phase': ['CIV_0', 'OVI_1']
    },
    use_lines=['MgII_2796', 'CII_1334', 'CIV_1548', 'OVI_1031']
)
```


### 2. Model Comparison

Compare different phase structures objectively:

```python
# Simple model
simple = cmbm.CMBMFitter(phases={'Single': ['CII_0']})
simple.load_data()
results_simple = simple.run_fit(output_dir='simple_model')

# Complex model
complex = cmbm.CMBMFitter(phases={'Low': ['SiII_0'], 'High': ['CIV_1']})
complex.load_data()
results_complex = complex.run_fit(output_dir='complex_model')

# Compare Bayesian evidence
bayes_factor = results_complex['logz'] - results_simple['logz']
print(f"Bayes factor: {bayes_factor:.2f}")
```

## Workflow Overview

```
┌─────────────────┐
│  Prepare Data   │
│  - VoigtFit     │
│  - Masks        │
│  - Grids        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Configure CMBM  │
│  - Define       │
│    phases       │
│  - Select lines │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Load Data     │
│  fitter.load    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Run Fit       │
│  Nested         │
│  Sampling       │
│  (1-24 hours)   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Analyze Results │
│  - Parameters   │
│  - Evidence     │
│  - Diagnostics  │
└─────────────────┘
```

## Core Concepts

### Phase Structure

A **phase** represents a distinct ionization state or velocity component:

```yaml
phases:
  Low_ionization:    # Phase name
    - "SiII_0"       # Cloud component (ion_id)
    - "MgII_1"
  High_ionization:
    - "CIV_0"
    - "SiIV_1"
```

Each component gets 5 fitted parameters (Z, nH, NHI, b_turb, z).

### Nested Sampling

CMBM uses UltraNest for:
- **Evidence calculation**: For model comparison
- **Posterior sampling**: Full parameter uncertainties
- **Multi-modal exploration**: Handles complex parameter spaces
- **Efficient sampling**: Focused on high-likelihood regions

### Velocity Masking

Exclude contaminated or irrelevant spectral regions:

```python
masks = {
    'CII_1334': [[-200, -50], [50, 200]],  
    'SiII_1190': [[-150, 150]]             
}
```

## Input Data Requirements

CMBM requires five input files:

1. **VoigtFit Dataset** (`.hdf5`)
   - Observed spectra with wavelength, flux, error
   - Line definitions and system redshift
   
2. **VoigtFit Results** (`.pickle`)
   - Initial parameter estimates
   - Used to set prior bounds
   
3. **Velocity Masks** (`.pkl`)
   - Define fitting regions
   - Exclude blends and contamination
   
4. **Cloudy Grids**
   - Precomputed photoionization models
   - Thin and thick grid sets
   
5. **Atomic Data** (`.dat`)
   - Transition wavelengths
   - Oscillator strengths
   - Damping parameters

See [Quick Start](quickstart.md) for details on preparing each file.

## Output

CMBM produces comprehensive output:

### Results Directory Structure

```
output_directory/
├── chains/
│   ├── equal_weighted_post.txt    # Posterior samples
│   └── weighted_post.txt          # Evidence-weighted samples
├── plots/
│   ├── run.pdf                    # Convergence diagnostics
│   ├── corner.pdf                 # Parameter correlations
│   └── trace.pdf                  # Parameter traces
├── results/
│   ├── params.json                # Parameter summaries
│   └── info.json                  # Run information
└── run.log                        # Detailed fitting log
```

### Key Outputs

**Parameter Posteriors**: Full uncertainty distributions
- Median, mean, standard deviation
- 16th/84th percentiles (1σ)
- Covariances and correlations

**Model Evidence**: log(Z) for Bayesian model comparison
- Higher evidence = better model fit
- Bayes factors for model selection

**Best-Fit Model**: Maximum likelihood parameters
- Model profiles for all fitted lines
- Residuals and χ² statistics

## Performance

### Typical Run Times

| System Complexity | Parameters | Lines | Cores | Approximate Time |
|-------------------|------------|-------|-------|------|
| Simple (2-3 phases) | 15 | 5-10 | 4 | 30 min |
| Moderate (4-5 phases) | 20 | 10-15 | 8 | 1 hr |
| Complex (7-10 phases) | 50 | 15-20 | 16 | 3 hrs |

### Optimization Tips

**Faster fitting:**
- Start with fewer lines (5-10) for testing
- Use fewer live points (200-400) initially
- Employ MPI parallelization
- Use simpler phase structures first

**Higher precision:**
- Increase live points (800-2000)
- Decrease dlogz target (0.1-0.3)
- Increase effective sample size (800-1000)
- More sampling steps (nsteps_factor = 3-4)

## Best Practices

### 1. Model Building Strategy

Start simple, add complexity gradually:

1. **Single phase** - Establish baseline
2. **Two phases** - Add ionization distinction
3. **Multiple components** - Add velocity structure
4. **Model comparison** - Justify complexity

### 2. Validation Checklist

Before trusting results:
- ✓ Convergence diagnostics look good
- ✓ No parameters hitting prior bounds
- ✓ Posteriors are unimodal (or multimodality is expected)
- ✓ Best-fit model visually matches data
- ✓ Parameters are physically reasonable
- ✓ Evidence is stable across runs

### 3. Workflow

1. **Exploratory analysis**: Fast, low-precision runs
2. **Model comparison**: Test different phase structures
3. **High-precision fit**: Final run with tight convergence
4. **Validation**: Check all diagnostics
5. **Documentation**: Save configurations and outputs

## Common Pitfalls

### ❌ Overfitting

**Problem**: Too many phases for available data

**Solution**: 
- Use Bayesian evidence for model comparison
- Start with simplest adequate model

### ❌ Poor Initial Values

**Problem**: VoigtFit results are far from true values

**Solution**:
- Manually adjust VoigtFit initial guesses
- Check VoigtFit convergence
- Use wider prior ranges if needed

### ❌ Insufficient Convergence

**Problem**: Evidence estimate not stable

**Solution**:
- Increase live points (min_num_live_points)
- Decrease dlogz target
- Let fit run longer
- Check diagnostics such as corner plots

### ❌ Velocity Masking Errors

**Problem**: Important absorption excluded or contamination included

**Solution**:
- Visually inspect all masked regions
- Be conservative (include slightly more than needed)
- Test sensitivity to mask boundaries


### MPI Parallelization

Speed up computation:

```bash
mpiexec -n 16 python fit_script.py
```

### Custom Grid Sets

Use different Cloudy grids:
- Different metallicity ranges
- Alternative ionizing spectra
- Custom density/column grids

## API Reference

For complete API documentation, see [api.md](api.md).

### Main Classes

**CMBMFitter**: Main fitting interface
- `load_data()`: Load all input files
- `run_fit()`: Execute nested sampling
- `get_results()`: Retrieve fitted parameters
- `summary()`: Print fitting summary

### Configuration

**Config**: Configuration management
- Load from YAML
- Override with parameters
- Save configurations

## Support and Contributing

### Getting Help

- **Documentation**: You're reading it!
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: Questions and community support
- **Email**: sameer@ou.edu

### Contributing

We welcome contributions! Areas of interest:
- Documentation improvements
- Bug fixes
- New features
- Example notebooks
- Performance optimizations

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## Citation

If you use CMBM in your research, please cite:

```bibtex
@software{cmbm2024,
  title={CMBM: Cloud-by-cloud Multiphase Bayesian ionization Modeling},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/cmbm}
}
```

## License

CMBM is released under the MIT License. See [LICENSE](../LICENSE) for details.

## Acknowledgments

CMBM builds upon:
- **VoigtFit** by J. K. Krogager
- **UltraNest** by Johannes Buchner
- **Cloudy** by Gary Ferland and collaborators

---

**Ready to get started?** Head to the [Installation Guide](installation.md) or [Quick Start Tutorial](quickstart.md)!

*Last updated: 2025*
