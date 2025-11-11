# CMBM: Cloud-by-cloud Multiphase Bayesian ionization Modeling

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**CMBM** is a Python package for fitting astronomical absorption line systems using nested sampling and Cloudy photoionization models. It enables robust Bayesian parameter estimation for complex multi-phase circumgalactic and intergalactic medium absorption systems.

## Table of Contents

- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Scientific Method](#scientific-method)
- [Configuration Guide](#configuration-guide)
- [Input Data](#input-data)
- [Running Fits](#running-fits)
- [Output](#output)
- [Post-processing and Analysis](#post-processing-and-analysis)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Key Features

- **Multi-phase modeling**: Fit absorption systems with multiple ionization phases simultaneously
- **Bayesian framework**: Full posterior sampling with UltraNest nested sampling
- **Cloudy integration**: Uses precomputed Cloudy photoionization model grids
- **YAML-based configuration**: Easy, reproducible setup via configuration files
- **Comprehensive post-processing**: Automated analysis tools for parameter extraction and visualization
- **HPC-ready**: Designed for parallel execution on high-performance computing systems

## Installation

```bash
pip install cmbm
```

Or install from source:
```bash
git clone https://github.com/sameeresque/cmbm.git
cd cmbm
pip install -e .
```
### Dependencies

CMBM will automatically install required packages:
- numpy, scipy, pandas, astropy
- ultranest (nested sampling)
- VoigtFit (spectral fitting)
- mpi4py (parallel computing)
- pyyaml (configuration)
- matplotlib (plotting)



## Quick Start

### 1. Create Configuration File

Copy and modify the example configuration:
```bash
cp config/default_config.yaml my_config.yaml
```

### Configuration File Example

`my_config.yaml`:

```yaml
# Default configuration for CMBM
# Cloud-by-cloud Multiphase Bayesian ionization Modeling

# ==============================================================================
# DATA PATHS - Update these paths for your data
# ==============================================================================
dataset_path: "path/to/your/dataset.hdf5"                          # Path to VoigtFit dataset (.hdf5)
voigtfit_path: "path/to/your/bestfit.pickle"                         # Path to VoigtFit results (.pickle)  
mask_path: "path/to/your/masks.pkl"                             # Path to velocity masks (.pkl)
transition_library_path: "path/to/atomdata.dat"              # Path to atomic transition data file

# Grid paths for Cloudy models
grid_paths:
  gridpath: "/path/to/cloudy/grids/"       # Thin grid directory
  gridpath_thick: "/path/to/cloudy/grids/" # Thick grid directory

grid_z: "0.2"                              # Grid redshift identifier

# ==============================================================================
# FITTING CONFIGURATION
# ==============================================================================

# Phase structure - defines the physical components to fit
# Each phase can contain multiple clouds with different ionization states
phases:
  Ph1: ["SiIII_0", "HI_1"]


# Spectral lines to include in the fit
use_lines: [
  "CII_1334", "CII_1036", "CIII_977", "FeII_1144", 
  "HI_1215", "HI_1025", "HI_972", "HI_949", "HI_937", "HI_930", "HI_926",
  "MgI_2852", "MgII_2796", "NII_1083", "NV_1238", "OI_1302", "OVI_1031",
  "SII_1253", "SIII_1012", "SiII_1190", "SiIII_1206", "SiIV_1393"
]

# Available ions in Cloudy grids (used for column density interpolation)
cloudy_ions: [
  "MgX", "AlII", "AlIII", "SIV", "CII", "HI", "CIV", "CIII", "SII", "TiII",
  "MgII", "NI", "OIII", "NaI", "SiII", "SiIII", "PIII", "NeVIII", "CaI", "FeI",
  "ZnII", "MnII", "SiIV", "SiI", "PI", "MgI", "NV", "CI", "SIII", "SVI",
  "FeII", "OIV", "CaII", "NII", "OII", "OI", "NIII", "AlI", "SI", "PII", "OV", "OVI"
]
# ==============================================================================
# MASKING CONFIGURATION
# ==============================================================================

masks:
  velocity_limits:
    lowlim: -160
    uplim: 190
    maxvel: 500
  custom_masks:
    CIII_977: [[-500, 5], [190, 500]]

# ==============================================================================
# NESTED SAMPLING PARAMETERS
# ==============================================================================

# Basic sampling settings
min_num_live_points: 400                    # Number of live points
dlogz: 0.5                                  # Target evidence accuracy (will be adjusted by ndim)
min_ess: 400                                # Minimum effective sample size
nsteps_factor: 2                            # Step sampler steps = nsteps_factor * ndim

# Advanced sampling settings
update_interval_volume_fraction: 0.2        # Volume update interval
max_num_improvement_loops: -1               # Maximum improvement loops (-1 = no limit)

# ==============================================================================
# OUTPUT SETTINGS
# ==============================================================================

output_dir: "path/to/your/SiIII0HI1"                  # Output directory name
resume: true                                # Resume previous run if available

# ==============================================================================
# GALAXY REDSHIFT
# ==============================================================================

# System redshift from VoigtFit dataset (z_sys)
# This is automatically loaded from the VoigtFit file

# Galaxy redshift for velocity reference in plots
# If not specified (null), z_sys will be used
zgal: 0.2284
# Example: zgal: 0.2284


# ==============================================================================
# PLOTTING CONFIGURATION
# ==============================================================================

plotting:
  # Velocity range for plots [vmin, vmax] in km/s
  velocity_range: [-30, 390]
  # Lines to exclude from model comparison plot
  deactivate_lines: ['CI_945','CI_1328','NI_1199','SIV_1062',
 'CaII_3969',
 'FeIII_1122',
 'FeII_1063',
 'FeII_1096',
 'FeII_1143',
 'FeII_2586',
 'NI_1201',
 'OI_1039',
 'SII_1259',
 'SI_1295',
 'SVI_944',
 'SiII_1020']

  
  # Lines to show as fully masked (grayed out)
  masked_lines: []



# ==============================================================================
# SLURM JOB SETTINGS
# ==============================================================================

slurm:
  job_name: "SiIII0HI1"           # Job name (-J flag)
  partition: "sooner_test"       # Partition name (--partition)
  container: "el9hw"             # Container image (--container)
  nodes: 1                       # Number of nodes (--nodes)
  ntasks_per_node: 8             # MPI tasks per node (--ntasks-per-node)
  mem_per_cpu: "16G"             # Memory per CPU (--mem-per-cpu)
  time: "04:00:00"               # Wall time (--time, HH:MM:SS)
```

### 2. Run the Fit
```python
import cmbm

# Load configuration and run fit
fitter = cmbm.CMBMFitter('my_config.yaml')
fitter.load_data()
results = fitter.run_fit()

# Print summary
fitter.summary()
```

### 3. Analyze Results

See [Post-processing and Analysis](#post-processing-and-analysis) section below.


## Scientific Method

CMBM implements a **cloud-by-cloud multiphase approach** to modeling absorption line systems:

1. **Phase Definition**: Define distinct ionization phases
2. **Cloud Components**: Each phase can contain multiple velocity components
3. **Physical Parameters**: Each cloud has 5 physical parameters (for the assumption of photoionization equilibrium and solar abundance pattern):
   - **Z**: Metallicity (relative to solar)
   - **nH**: Hydrogen density (cm⁻³)
   - **NHI**: Neutral hydrogen column density (cm⁻²)
   - **b_turb**: Turbulent broadening (km/s)
   - **z**: Redshift

4. **Cloudy Models**: Column densities computed from precomputed Cloudy photoionization grids
5. **Bayesian Fitting**: UltraNest nested sampling provides full posterior distributions
6. **Model Comparison**: Bayesian evidence enables robust model comparison

### Derived Quantities

From the 5 fitted parameters, CMBM computes:
- **Temperature** (T): From Cloudy models run with the assumption of photoionization thermal equilibrium
- **Cloud size** (L): From total hydrogen column density and gas volume density
- **Total hydrogen column** (NH): From Cloudy models
- **Velocity**: Relative to galaxy redshift
- **Ion-specific b-parameters**: Including thermal broadening

### Bayesian Framework

CMBM uses a Bayesian approach to fit absorption line systems:

#### Setup

The observed data **D** = {d_i} are measurements of an underlying model function f(λ; **θ**) at N spectral pixels, where **θ** = {Z, n_H, N(HI), b_nt, z} are the physical parameters for each cloud in photoionization equilibrium.

The relationship between data and model:
```
d_i = f_i(θ) + ε_i,    ε_i ~ N(0, σ_i²)
```

where ε_i represents Gaussian measurement uncertainties with standard deviation σ_i.

#### Likelihood Function

The likelihood of observing the data given the model parameters is:
```
p(D|θ, M) = ∏[i=1 to N] (1/(σ_i√(2π))) exp[-(1/2)((d_i - f_i(θ))/σ_i)²]

           ∝ exp[-(1/2) Σ_i ((d_i - f_i(θ))/σ_i)²]
           
           = exp[-χ²(θ)/2]
```

where χ²(θ) is the chi-squared statistic measuring goodness of fit.

#### Posterior Distribution

Using Bayes' theorem with prior density π(θ):
```
p(θ|D, M) ∝ π(θ) × p(D|θ, M)
          ∝ π(θ) exp[-χ²(θ)/2]
```

The posterior distribution combines prior knowledge with observational constraints.

#### Parameter Inference

This is accomplished using UltraNest, which efficiently:
- Handle multimodal posterior distributions
- Compute the Bayesian evidence (log Z)
- Generate posterior samples for parameter estimation
- Provide robust uncertainty quantification

#### Advantages of Bayesian Approach

- **Full uncertainty propagation**: Posterior distributions capture parameter degeneracies and correlations
- **Model comparison**: Bayesian evidence enables rigorous comparison of different phase structures
- **Physical constraints**: Prior distributions incorporate known physical bounds (e.g., density ranges, metallicity limits)
- **Multimodal solutions**: Nested sampling naturally handles multiple parameter solutions

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
   - Contains wavelength, flux, error arrays for each absorption line
   - Created using VoigtFit's `dataset.save()` method

2. **VoigtFit Results** (`.pickle`): Initial parameter estimates from VoigtFit
   - Provides starting guesses for cloud redshifts and velocities
   - Created by VoigtFit's fitting procedure

3. **Velocity Masks** (`.pkl`): Velocity regions to include/exclude in fits
   - Dictionary specifying masked velocity ranges per line
   - Can be created with custom scripts or VoigtFit utilities

4. **Cloudy Grids**: Precomputed photoionization model grids
   - Must match `grid_z` parameter in config

5. **Atomic Data** (`.dat`): Transition wavelengths, oscillator strengths, damping constants
   - Standard format atomic line database

## Running Fits

### Interactive Python Session
```python
import cmbm

fitter = cmbm.CMBMFitter('my_config.yaml')
fitter.load_data()
results = fitter.run_fit()
```

### Command Line
```bash
python -c "import cmbm; f = cmbm.CMBMFitter('my_config.yaml'); f.load_data(); f.run_fit()"
```

### SLURM Batch Script
```bash
#!/bin/bash
#SBATCH --job-name=cmbm_fit
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

module load Python/3.11.3-GCCcore-12.3.0
module load OpenMPI/4.1.4-GCC-12.2.0

cd /path/to/your/work/directory

python << EOF
import cmbm
fitter = cmbm.CMBMFitter('my_config.yaml')
fitter.load_data()
fitter.run_fit()
EOF
```

### Parallel Execution with MPI
```bash
mpiexec -n 8 python -c "import cmbm; f = cmbm.CMBMFitter('my_config.yaml'); f.load_data(); f.run_fit()"
```

## Output

CMBM produces comprehensive output including:

- **Nested sampling chains**: Full posterior samples
- **Parameter summaries**: Medians, uncertainties, correlations
- **Model comparison**: Bayesian evidence for different phase structures
- **Best-fit models**: Model profiles at maximum likelihood parameters
- **Diagnostic plots**: Residuals, corner plots, trace plots

### Directory Structure
```
output_dir/
├── chains/
│   ├── equal_weighted_post.txt    # Equal-weighted posterior samples
│   ├── weighted_post.txt          # Weighted posterior samples
│   └── ...
├── info/
│   ├── results.json               # Summary statistics
│   ├── post_summary.json          # Parameter summaries
│   └── ...
└── plots/                         # Diagnostic plots (if enabled)
```

### Key Output Files

- **`equal_weighted_post.txt`**: Posterior samples with equal weights (for analysis)
- **`weighted_post.txt`**: Posterior samples with importance weights
- **`results.json`**: Contains:
  - Log evidence (log Z)
  - Evidence uncertainty
  - Number of likelihood evaluations
  - Sampling efficiency metrics
- **`post_summary.json`**: Parameter statistics (median, mean, std, etc.)

## Post-processing and Analysis

After fitting, analyze your results with the post-processing tool:

### Basic Usage
```bash
python postprocess/analyze_results.py \
    --config my_config.yaml \
    --output_dir my_results \
    --save_dir analysis_output \
    --make_plots
```

### Command-Line Options

#### Required Arguments
- `--config`: Path to configuration file used for fitting
- `--output_dir`: Directory containing CMBM output (chains/, info/, etc.)
- `--save_dir`: Directory where analysis results will be saved

#### Optional Flags
- `--make_plots`: Generate model comparison plots
- `--skip_derived`: Skip computing derived quantities (load from file instead)
- `--skip_stats`: Skip computing statistics (load from file instead)
- `--skip_corner`: Skip corner plot generation
- `--n_processes`: Number of parallel processes (default: 12)

### Skip Flags Explained

The skip flags allow you to avoid recomputing time-consuming steps when re-running analysis:

**`--skip_derived`**: Skips computation of derived quantities (temperature, cloud size, velocities, ion columns)
- Use when: You have already run the analysis once and only want to update plots
- Requires: `posterior_full.pkl` must exist in `save_dir`

**`--skip_stats`**: Skips computation of parameter statistics (medians, uncertainties)
- Use when: You've already computed statistics and only need plots
- Requires: `parameter_statistics.pkl` and `ion_properties.pkl` must exist

**`--skip_corner`**: Skips generating corner plots
- Use when: You only want the model comparison plot, not parameter correlations

### Model Comparison Plot

The model comparison plot (`model_comparison.pdf`) shows:

**Left panels**: Absorption line profiles with:
- Gray points: Observed data with error bars
- Colored lines: Individual cloud models
- Black line: Combined model (all clouds)
- Gray shading: Masked regions
- Vertical colored lines: Cloud velocities

**Right panel**: Parameter distributions as violin plots showing:
- Filled violins: 1σ posterior distributions
- Dashed outlines: 3σ posterior distributions
- Black stars: Median values
- Parameters shown: Z, nH, NHI, T, L, b_turb

Lines are automatically ordered:
1. HI lines first, sorted by f*λ (oscillator strength × wavelength)
2. Metal lines alphabetically by element, then by ionization state


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

We welcome contributions! Please contact the author.

## License

CMBM is released under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

CMBM builds upon:
- [VoigtFit](https://github.com/jkrogager/VoigtFit) for spectral line fitting
- [UltraNest](https://github.com/JohannesBuchner/UltraNest) for nested sampling
- [Cloudy](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home) photoionization code

---

**Happy fitting!** 🌟
