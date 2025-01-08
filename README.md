# ASE Cantera Microkinetics

**ASE Cantera Microkinetics** is a Python library that facilitates microkinetic modeling in heterogeneous catalysis. It bridges the gap between atomistic simulations using [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) and reactor simulations via [Cantera](https://cantera.org/). By integrating these tools, the library enables streamlined workflows for studying reaction mechanisms, rate equations, and catalytic performance under realistic operating conditions.

## Features

- **Reaction Mechanism Definition**: Easily define elementary reaction steps with support for custom reaction networks.
- **Thermochemistry Data Integration**: Seamlessly incorporate thermodynamic properties from ASE-calculated data or other sources.
- **Dynamic Simulations**: Perform transient and steady-state reactor simulations using Cantera's robust solvers.
- **Customizable Workflows**: Extendable and modular architecture for specific catalysis problems.

## Installation

To install the package, clone the repository and use the following commands:

```bash
git clone https://github.com/raffaelecheula/ase_cantera_microkinetics.git
cd ase_cantera_microkinetics
pip install -e .
```

## Usage

An example of how to use **ase_cantera_microkinetics** is in `ase_cantera_microkinetics\examples\microkinetics_integration.py`

## Requirements

- Python 3.8 or later
- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [Cantera](https://cantera.org/)
- NumPy, SciPy, and Matplotlib (for data analysis and visualization)

## Contributing

Contributions are welcome! Please fork the repository, create a feature branch, and submit a pull request. Ensure your code follows PEP 8 guidelines.

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Authors and Acknowledgments

Developed by Raffaele Cheula. Special thanks to the open-source community and the developers of ASE and Cantera for their foundational contributions.

## Support

Please open a GitHub Issue or contact the maintainer at [cheula.raffaele@gmail.com] for issues or questions.
