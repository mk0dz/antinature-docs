# Antinature

A Python module for quantum chemical calculations of antimatter and mixed matter-antimatter systems.

## Overview

Antinature extends conventional quantum chemistry methods to accurately model positrons, positronium, and complex antimatter structures, providing researchers with tools to explore antimatter chemistry and physics.

Key features:
- Specialized basis sets for positrons
- Extended Hamiltonians with annihilation operators
- Self-consistent field methods for antimatter systems
- Relativistic corrections for accurate antimatter modeling
- Tools for calculating annihilation rates and lifetimes

## Installation

```bash
# Install from PyPI
pip install antinature

# For development installation
git clone https://github.com/yourusername/antinature.git
cd antinature
pip install -e .
```

Prerequisites:
- Python 3.7+
- NumPy
- SciPy
- Matplotlib
- Qiskit (optional, for quantum computing functionality)

## Quick Start

```python
import numpy as np
from antinature.core import MolecularData, MixedMatterBasis
from antinature.core.scf import AntinatureSCF
from antinature.core.hamiltonian import AntinatureHamiltonian
from antinature.core.integral_engine import AntinatureIntegralEngine
from antinature.specialized import PositroniumSCF

# Create a positronium system
pos_system = MolecularData.positronium()

# Create a specialized basis set
basis = MixedMatterBasis()
basis.create_positronium_basis(quality='extended')

# Set up integral engine
integral_engine = AntinatureIntegralEngine(use_analytical=True)
basis.set_integral_engine(integral_engine)

# Build Hamiltonian
hamiltonian = AntinatureHamiltonian(
    molecular_data=pos_system,
    basis_set=basis,
    integral_engine=integral_engine,
    include_annihilation=True
)
hamiltonian_result = hamiltonian.build_hamiltonian()

# Run SCF calculation
scf = PositroniumSCF(
    hamiltonian=hamiltonian_result,
    basis_set=basis,
    molecular_data=pos_system,
    max_iterations=100,
    convergence_threshold=1e-6
)
scf.solve_scf()

# Get results
energy = scf.compute_energy()
print(f"Positronium ground state energy: {energy:.6f} Hartree")
```

## Documentation

Comprehensive documentation is available in the `docs/` directory and includes:

- [API Reference](docs/api_reference.md)
- [Theory Guide](docs/theory.md)
- [Examples](docs/examples.md)

## Tutorials

Step-by-step tutorials are available in the `Tutorials/` directory:

1. [Introduction to Antimatter Calculations](Tutorials/01_intro_to_antimatter.ipynb)
2. [Working with Positronium](Tutorials/02_working_with_positronium.ipynb)
3. [Advanced Basis Sets](Tutorials/03_advanced_basis_sets.ipynb)
4. [Relativistic Effects](Tutorials/04_relativistic_effects.ipynb)
5. [Quantum Computing Applications](Tutorials/05_quantum_computing.ipynb)

## Examples

Additional examples demonstrating specific applications are in the `examples/` directory:

- Simple antimatter atoms (He+, H+)
- Molecules with positrons
- Annihilation rate calculations
- Positronium excited states

## Citation

If you use Antinature in your research, please cite:

```
@article{antinature2023,
  title={Antinature: A Quantum Chemistry Framework for Antimatter Systems},
  author={Author, A. and Author, B.},
  journal={Journal of Computational Chemistry},
  year={2023},
  volume={},
  pages={}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The pioneering work of Paul Dirac on antimatter theory
- Experimental advances at CERN and other antimatter research facilities
- Contributors and supporters of the Antinature project
