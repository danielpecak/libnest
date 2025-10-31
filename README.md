# libnest

**Library for Neutron Star Physics**

A Python library for nuclear matter and neutron star physics calculations, implementing the Brussels-Montreal (BSk) energy density functional and related physical models.

[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## 🔬 Overview

`libnest` provides computational tools for:

- **Energy density functionals**: Brussels-Montreal Skyrme parametrization (BSk31)
- **Nuclear matter properties**: Equations of state, effective masses, pairing fields
- **Neutron star physics**: Inner crust structure, superfluid properties, pasta phases
- **Unit conversions**: Nuclear physics units (MeV, fm⁻³, etc.) and astronomical units
- **Visualization**: Ready-to-use plotting functions for physical quantities

This library is designed for researchers working on:
- Dense matter in neutron stars
- Nuclear superfluidity and pairing
- Neutron star crust composition
- Nuclear energy density functionals

## 📦 Installation

### From source (recommended for development)

```bash
git clone https://github.com/danielpecak/libnest.git
cd libnest
pip install -e .
```

### Requirements

- Python ≥ 3.7
- NumPy
- Matplotlib
- Pandas
- SciPy (optional, for advanced features)

See `requirements.txt` for specific versions.

## 🚀 Quick Start

```python
import libnest
import libnest.bsk
import libnest.definitions

# Example 1: Unit conversion
rho_fm3 = 0.16  # nuclear saturation density
rho_gcm3 = libnest.units.fm3togcm3(rho_fm3)
print(f"Density: {rho_fm3} fm⁻³ = {rho_gcm3:.2e} g/cm³")

# Example 2: Calculate Fermi energy
rho_n = 0.08  # neutron density in fm⁻³
kF = libnest.definitions.rho2kf(rho_n)
eF = libnest.definitions.eF_n(kF)
print(f"Fermi energy: {eF:.2f} MeV")

# Example 3: BSk energy functional
rho_n = rho_p = 0.08  # symmetric nuclear matter
E_per_A = libnest.bsk.energy_per_nucleon(rho_n, rho_p)
print(f"Energy per nucleon: {E_per_A:.2f} MeV")

# Example 4: Neutron pairing gap
delta_n = libnest.bsk.neutron_pairing_field(rho_n)
print(f"Pairing gap: {delta_n:.2f} MeV")
```

Run the included example:
```bash
python main.py
```

## 📚 Documentation

Full documentation is available at: [https://libnest.readthedocs.io](https://libnest.readthedocs.io) *(coming soon)*

### Building documentation locally

```bash
cd docs
make html
```

The documentation will be in `docs/_build/html/`.

## 🧪 Modules

### Core Modules

- **`units`**: Physical constants and unit conversion functions
- **`definitions`**: Basic definitions (Fermi wavevector, energy, velocities)
- **`bsk`**: Brussels-Montreal Skyrme functional (BSk31 parametrization)
- **`tools`**: Utility functions (center of mass, condensation energy, slicing)
- **`myio`**: Input/output functions for data files

### Analysis Modules

- **`plots`**: Plotting functions for uniform matter properties
- **`real_data_plots`**: Visualization for simulation output
- **`delta_and_temperature`**: Temperature-dependent pairing analysis
- **`pasta`**: Nuclear pasta phases in the crust
- **`nucleus`**: Nuclear structure calculations

## 📖 Examples

See the `examples/` directory for more comprehensive examples:

- `examples/tools_example.py` - Using utility functions
- `examples/legacy_tests.py` - Legacy test cases and plotting examples

## 🧮 Physical Constants

The library uses natural units with ℏ = c = 1 where convenient. Key constants:

| Constant | Value | Unit | Description |
|----------|-------|------|-------------|
| `HBARC` | 197.327 | MeV·fm | ℏc |
| `MN` | 939.565 | MeV | Neutron mass |
| `MP` | 938.272 | MeV | Proton mass |
| `RHOSAT` | 0.16 | fm⁻³ | Nuclear saturation density |

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

**Daniel Pęcak**
- Email: daniel.pecak@pw.edu.pl
- Affiliation: Warsaw University of Technology, Université Libre de Bruxelles
- On leave from: Institute of Physics, Polish Academy of Sciences

## 🙏 Acknowledgments

This work has been supported by:
- Polish National Science Centre (NCN)
- PLGrid Infrastructure
- LUMI Supercomputer (EuroHPC)

## 📚 References

If you use this library in your research, please cite:

- S. Goriely, N. Chamel, J.M. Pearson, *Phys. Rev. C* **88**, 024308 (2013) - BSk functional
- N. Chamel, S. Goriely, J.M. Pearson, *Phys. Rev. C* **80**, 065804 (2009) - Pairing model

See `docs/bibliography.rst` for a complete list of references.

## 🔗 Related Projects

- [SkyNET](https://github.com/nuclear-physics/skynet) - Nuclear reaction network
- [CompOSE](https://compose.obspm.fr/) - Equations of state database

---

**Version**: 2024.01  
**Last Updated**: November 2025