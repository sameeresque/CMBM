# CMBM: Cloud-by-cloud Multiphase Bayesian ionization Modeling

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**CMBM** is a Python package for fitting astronomical absorption line systems using nested sampling and Cloudy photoionization models. It enables robust Bayesian parameter estimation for complex multi-phase circumgalactic and intergalactic medium absorption systems.

## Key Features

- **Multi-phase modeling**: Fit absorption systems with multiple ionization phases simultaneously
- **Bayesian framework**: Full posterior sampling with UltraNest nested sampling
- **Cloudy integration**: Uses precomputed Cloudy photoionization model grids
- **Flexible configuration**: Easy setup via YAML configuration files or direct parameters

## Quick Start

### Installation

```bash
pip install cmbm
```

Or install from source:
```bash
git clone https://github.com/sameeresque/cmbm.git
cd cmbm
pip install -e .
```

### Basic Usage

```python
import cmbm

# Method 1: Using configuration file
fitter = cmbm.CMBMFitter('my_config.yaml')
fitter.load_data()
results = fitter.run_fit()

# Method 2: Direct parameter specification
fitter = cmbm.CMBMFitter(
    dataset_path='path/to/dataset.hdf5',
    voigtfit_path='path/to/bestfit.pickle',
    mask_path='path/to/masks.pkl',
    grid_paths={'gridpath': '/path/to/cloudy/grids/'},
    transition_library_path='path/to/atomdata.dat',
    phases={'Ph1': ['SiIII_0', 'HI_1']},
    use_lines=['CII_1334', 'HI_1215', 'MgII_2796']
)
fitter.load_data()
results = fitter.run_fit(output_dir='my_results')

# Get results
summary = fitter.get_results()
print(f"Log evidence: {summary['log_evidence']:.2f}")

# Print parameter summaries
fitter.summary()
```

### Configuration File Example

Create a configuration file `my_config.yaml`:

```yaml
# Data paths
dataset_path: "path/to/your/dataset.hdf5"
voigtfit_path: "path/to/your/bestfit.pickle"
mask_path: "path/to/your/masks.pkl"
transition_library_path: "path/to/atomdata.dat"

grid_paths:
  gridpath: "/path/to/cloudy/grids/"
  gridpath_thick: "/path/to/cloudy/grids/"
grid_z: "0.2"

# Define your absorption system phases
phases:
  Low_ionization: ["SiII_0", "MgII_1", "FeII_2"]
  High_ionization: ["CIV_0", "SiIV_1", "OVI_2"]

# Select lines to fit
use_lines:
  - "CII_1334"
  - "SiII_1190"
  - "MgII_2796"
  - "CIV_1548"
  - "SiIV_1393"
  - "HI_1215"

# Sampling settings
min_num_live_points: 800
output_dir: "my_absorption_system"
```

## Scientific Method

CMBM implements a **cloud-by-cloud multiphase approach** to modeling absorption line systems:

1. **Phase Definition**: Define distinct ionization phases
2. **Cloud Components**: Each phase can contain multiple velocity components
3. **Physical Parameters**: Each cloud has 5 physical parameters:
   - **Z**: Metallicity (relative to solar)
   - **nH**: Hydrogen density (cm⁻³)
   - **NHI**: Neutral hydrogen column density (cm⁻²)
   - **b_turb**: Turbulent broadening (km/s)
   - **z**: Redshift

4. **Cloudy Models**: Column densities computed from precomputed Cloudy photoionization grids
5. **Bayesian Fitting**: UltraNest nested sampling provides full posterior distributions
6. **Model Comparison**: Bayesian evidence enables robust model comparison

## Requirements

- Python ≥ 3.7
- numpy ≥ 1.24.4
- scipy ≥ 1.10.1
- pandas ≥ 2.0.3
- astropy ≥ 5.2.2
- ultranest ≥ 4.4.0
- VoigtFit ≥ 3.21.3
- mpi4py ≥ 3.1.5
- pyyaml ≥ 5.4.1

## Input Data

CMBM requires several input files:

1. **VoigtFit Dataset** (`.hdf5`): Observed spectra prepared with VoigtFit
2. **VoigtFit Results** (`.pickle`): Initial parameter estimates from VoigtFit
3. **Velocity Masks** (`.pkl`): Velocity regions to include/exclude in fits
4. **Cloudy Grids**: Precomputed photoionization model grids
5. **Atomic Data** (`.dat`): Transition wavelengths, oscillator strengths, etc.

## Output

CMBM produces comprehensive output including:

- **Nested sampling chains**: Full posterior samples
- **Parameter summaries**: Medians, uncertainties, correlations
- **Model comparison**: Bayesian evidence for different phase structures
- **Best-fit models**: Model profiles at maximum likelihood parameters
- **Diagnostic plots**: Residuals, corner plots, trace plots

## Examples

### Simple Two-Phase System

```python
# Model a simple low/high ionization system
fitter = cmbm.CMBMFitter(
    dataset_path='quasar_spectrum.hdf5',
    phases={
        'Low_ion': ['SiII_0', 'MgII_1'],  
        'High_ion': ['CIV_0', 'OVI_1']
    },
    use_lines=['SiII_1190', 'MgII_2796', 'CIV_1548', 'OVI_1031']
)
```

## Citation

If you use CMBM in your research, please consider citing:

```bibtex
@ARTICLE{2024MNRAS.530.3827S,
       author = {{Sameer} and {Charlton}, Jane C. and {Wakker}, Bart P. and {Kacprzak}, Glenn G. and {Nielsen}, Nikole M. and {Churchill}, Christopher W. and {Richter}, Philipp and {Muzahid}, Sowgat and {Ho}, Stephanie H. and {Nateghi}, Hasti and {Rosenwasser}, Benjamin and {Narayanan}, Anand and {Ganguly}, Rajib},
        title = "{Cloud-by-cloud multiphase investigation of the circumgalactic medium of low-redshift galaxies}",
      journal = {\mnras},
     keywords = {methods: statistical, galaxies: evolution, galaxies: haloes, galaxies: individual...intergalactic medium, quasars: absorption lines, Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = jun,
       volume = {530},
       number = {4},
        pages = {3827-3854},
          doi = {10.1093/mnras/stae962},
archivePrefix = {arXiv},
       eprint = {2403.05617},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.3827S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
@ARTICLE{2021MNRAS.501.2112S,
       author = {{Sameer} and {Charlton}, Jane C. and {Norris}, Jackson M. and {Gebhardt}, Matthew and {Churchill}, Christopher W. and {Kacprzak}, Glenn G. and {Muzahid}, Sowgat and {Narayanan}, Anand and {Nielsen}, Nikole M. and {Richter}, Philipp and {Wakker}, Bart P.},
        title = "{Cloud-by-cloud, multiphase, Bayesian modelling: application to four weak, low-ionization absorbers}",
      journal = {\mnras},
     keywords = {quasars: general, quasars: absorption lines, galaxies: general, galaxies: evolution, Astrophysics - Astrophysics of Galaxies},
         year = 2021,
        month = feb,
       volume = {501},
       number = {2},
        pages = {2112-2139},
          doi = {10.1093/mnras/staa3754},
archivePrefix = {arXiv},
       eprint = {2012.00021},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.2112S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

CMBM is released under the MIT License. See [LICENSE](LICENSE) for details.

## Support

- **Documentation**: [https://cmbm.readthedocs.io/](https://cmbm.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/yourusername/cmbm/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/cmbm/discussions)

## Acknowledgments

CMBM builds upon:
- [VoigtFit](https://github.com/jkrogager/VoigtFit) for spectral line fitting
- [UltraNest](https://github.com/JohannesBuchner/UltraNest) for nested sampling
- [Cloudy](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home) photoionization code

---

**Happy fitting!** 🌟
