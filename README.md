# StereoMolGraph

[![View Documentation](https://img.shields.io/badge/📖-Documentation-8CA1AF)](https://stereomolgraph.readthedocs.io)
[![Python Versions](https://img.shields.io/badge/Python-3.10|3.11|3.12|3.13-blue?logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?logo=opensourceinitiative)](https://opensource.org/licenses/MIT)
![PyPI - Version](https://img.shields.io/pypi/v/StereoMolGraph?link=https%3A%2F%2Fpypi.org%2Fproject%2FStereoMolGraph)
[![Documentation Status](https://readthedocs.org/projects/stereomolgraph/badge/?version=dev)](https://stereomolgraph.readthedocs.io/en/latest/?badge=dev)
[![GitHub](https://img.shields.io/badge/GitHub-View%20on%20GitHub-blue?logo=github)](https://github.com/maxim-papusha/StereoMolGraph)
[![Type Check](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_type_check.yaml/badge.svg?branch=main&event=push)](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_type_check.yaml)
[![Test Examples](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_example_test.yaml/badge.svg?branch=main&event=push)](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_example_test.yaml)
[![Unit Tests](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_unit_test.yaml/badge.svg?branch=main&event=push)](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_unit_test.yaml)
[![DOI](https://zenodo.org/badge/999259345.svg)](https://doi.org/10.5281/zenodo.16360310)
[![PyPI Downloads](https://static.pepy.tech/badge/stereomolgraph)](https://pepy.tech/projects/stereomolgraph)

StereoMolGraph (SMG) is a lightweight Python library for representing molecules and transition states with explicit, local stereochemistry. It provides:

- Graph types for molecules and reactions (with/without stereo and stereo changes)
- Includes non tetrahedral stereocenters and changing stereochemistry in reactions
- Fast approximate hashing via Weisfeiler–Lehman color refinement
- Robust equality/isomorphism via a VF2++-style algorithm extended for stereochemistry and reactions
- Bidirectional conversion from / to RDKit
- Construction from 3D coordinates with automatic local stereo inference


## Design philosophy

- Unopinionated about bond orders, charge and electronic state
- SMG focuses on the connectivity and stereochemistry. 
- Stereochemistry describes relative spatial arrangement. No absolute stereochemistry.
- Transparent: Simple 2D visualization in IPython notebooks


## Examples

Explore the working example notebooks in `docs/examples/` (executed in CI). For rendered examples and guides, see the documentation: https://stereomolgraph.readthedocs.io

## RDKit interoperability notes

- Hydrogens must be explicit for stereo-safe bidirectional conversion.
- Supports tetrahedral and non tetrahedral stereochemistry during conversion.
- Bond orders, charges, unpaired electrons and other properties are not used!

## Installation

Install from PyPI:

```bash
pip install stereomolgraph
```

## Feedback and support

Bug reports and feature requests are welcome — please open an issue on GitHub:

- Issues: https://github.com/maxim-papusha/StereoMolGraph/issues

If you have questions or ideas that don’t fit a template, you can still open an issue and tag it appropriately.

## Citation

If you use StereoMolGraph in your work, please cite the Zenodo record:

[![DOI](https://zenodo.org/badge/999259345.svg)](https://doi.org/10.5281/zenodo.16360310)

## License

MIT License — see `LICENSE`.

