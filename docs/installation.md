# Installation Guide

This guide covers installing CMBM and all its dependencies on various platforms.

## Quick Installation

For most users, installation is straightforward:

```bash
pip install cmbm
```

That's it! This will install CMBM and all required dependencies.

## System Requirements

### Minimum Requirements

- **Python**: 3.7 or higher
- **Operating System**: Linux, macOS, or Windows (Linux recommended for HPC)
- **Memory**: 8 GB RAM minimum (16+ GB recommended for complex systems)
- **Storage**: 10 GB for grids + output (depends on grid size)
- **CPU**: Multi-core processor (8+ cores recommended for nested sampling with MPI)


## Detailed Installation Instructions

### Step 1: Set Up Python Environment

We strongly recommend using a virtual environment or conda environment:

#### Option A: Using venv (built-in)

```bash
# Create virtual environment
python3 -m venv cmbm_env

# Activate environment (Linux/macOS)
source cmbm_env/bin/activate

# Activate environment (Windows)
cmbm_env\Scripts\activate
```

#### Option B: Using conda

```bash
# Create conda environment
conda create -n cmbm python=3.10

# Activate environment
conda activate cmbm
```

### Step 2: Install CMBM

#### From PyPI (Recommended)

```bash
pip install cmbm
```

#### From Source (Development)

```bash
# Clone repository
git clone https://github.com/yourusername/cmbm.git
cd cmbm

# Install in development mode
pip install -e .

# Or install with development dependencies
pip install -e .[dev,docs]
```

### Step 3: Verify Installation

```python
# Test import
import cmbm
print(f"CMBM version: {cmbm.__version__}")

# Check dependencies
import numpy
import scipy
import pandas
import astropy
import ultranest
import VoigtFit
import mpi4py

print("All dependencies installed successfully!")
```

## Dependencies

CMBM has the following core dependencies:

### Required Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥1.24.4 | Array operations |
| scipy | ≥1.10.1 | Scientific computing |
| pandas | ≥2.0.3 | Data structures |
| astropy | ≥5.2.2 | Astronomical calculations |
| ultranest | ≥4.4.0 | Nested sampling |
| VoigtFit | ≥3.21.3 | Voigt profile fitting |
| mpi4py | ≥3.1.5 | MPI parallelization |
| pyyaml | ≥5.4.1 | Configuration files |
| roman | ≥3.0 | Roman numeral conversion |

### Grid Directory Structure

Organize grids as:

```
/path/to/grids/
├── df_piethin/
│   ├── dfpiethin_0.2.pkl
│   ├── dfpiethin_0.2.pkl_coldens_HI
│   ├── dfpiethin_0.2.pkl_coldens_CII
│   ├── dfpiethin_0.2.pkl_coldens_CIV
│   └── ... (one file per ion + logT + logNHtot)
└── df_piethick/
    ├── df_pie_thick_0.2.pkl
    ├── df_pie_thick_0.2.pkl_coldens_HI
    └── ... (similar structure)
```

### Grid Size

Typical grid set:
- **Size**: 2-10 GB per redshift
- **Ions**: 40+ species
- **Parameters**: Metallicity, density, column density, temperature

## Testing Your Installation

### Basic Test

```python
import cmbm

# Create a fitter (will use defaults)
fitter = cmbm.CMBMFitter()
print("CMBM initialized successfully!")
```

### Comprehensive Test

```python
# Test all major imports
import numpy as np
import scipy
import pandas as pd
from astropy import constants
import ultranest
import VoigtFit
from mpi4py import MPI
import yaml

print("✓ All dependencies imported successfully")

# Test CMBM functionality
import cmbm
fitter = cmbm.CMBMFitter()
print("✓ CMBM fitter created successfully")

# Test MPI
comm = MPI.COMM_WORLD
print(f"✓ MPI working: {comm.size} process(es)")

print("\nInstallation test complete! ✓")
```

## Getting Help

If you encounter installation issues:

1. Check this documentation thoroughly
2. Search [GitHub Issues](https://github.com/sameeresque/cmbm/issues)
3. Ask on [GitHub Discussions](https://github.com/sameeresque/cmbm/discussions)
4. Contact developers: contact@yourdomain.com

## Next Steps

After successful installation:

1. Read the [Quick Start Guide](quickstart.md)
2. Try the [Basic Examples](../examples/basic_usage.py)
3. Configure your [data paths](quickstart.md#preparing-your-data)

