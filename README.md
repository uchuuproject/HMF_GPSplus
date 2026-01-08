# HMF: Halo Mass Function Model (Fernandez-Garcia 2025)

This repository contains the implementation of the Halo Mass Function (HMF) model based on the Fernandez-Garcia et al. (2025) formalism. The code is designed to work with the `Colossus` cosmology toolkit and provides a robust way to calculate halo abundances across different mass ranges and cosmologies.

## Installation

### Prerequisites

You will need Python 3.7+ and the following dependencies:

- `numpy`
- `scipy`
- `colossus`

### Install Dependencies

```bash
pip install numpy scipy colossus
```

## Usage

You can use the `HaloMassFunction` class to calculate the mass function at a specific redshift.

```python
import numpy as np
from HMF import HaloMassFunction

# Initialize the model at z=0 with a specific mass definition
hmf_model = HaloMassFunction(z=0.0, mdef='m200b')

# Define a range of masses (M_sun/h)
masses = np.logspace(10, 15, 100)

# Calculate dn/dlnM
abundance = hmf_model.n0(masses)

# abundance is now an array containing the halo density
```

## Model Details

The model implements the following components:
- **b(M) and c(M)**: Mass-dependent fitting functions for the peak height logic.
- **F(M)**: The cumulative collapse fraction.
- **n0(M)**: The differential mass function $dn/d\ln M$.

## References

If you use this code in your research, please cite:
*[Fernandez-Garcia+25](https://arxiv.org/abs/2512.05847)*.

## License

This project is licensed under the MIT License.
