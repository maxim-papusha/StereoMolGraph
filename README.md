<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/_static/img/logo_smg_dark.svg">
  <img alt="Logo" src="docs/_static/img/logo_smg.svg">
</picture>

# StereoMolGraph #

[![PyPI](https://img.shields.io/pypi/v/StereoMolGraph?style=flat-square&logo=pypi&logoColor=white&color=3775A9)](https://pypi.org/project/StereoMolGraph/)
[![Python](https://img.shields.io/pypi/pyversions/StereoMolGraph?style=flat-square&logo=python&logoColor=white&color=3776AB&label=Python)](https://pypi.org/project/StereoMolGraph/)
[![License: MIT](https://img.shields.io/badge/License-MIT-d4a017?style=flat-square&logo=opensourceinitiative&logoColor=white)](https://opensource.org/licenses/MIT)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/maxim-papusha/StereoMolGraph/run_unit_test.yaml?branch=main&style=flat-square&label=tests)](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_unit_test.yaml)

[![Documentation](https://img.shields.io/badge/Documentation-docs-4C6A92?style=flat-square&logo=readthedocs&logoColor=white)](https://stereomolgraph.readthedocs.io)
[![DOI](https://img.shields.io/badge/DOI-10.26434%2Fchemrxiv--2025--0g4wn-ffcc00?style=flat-square)](https://chemrxiv.org/doi/full/10.26434/chemrxiv-2025-0g4wn)

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

[![DOI](https://img.shields.io/badge/DOI-10.26434%2Fchemrxiv--2025--0g4wn-ffcc00?style=flat-square)](https://chemrxiv.org/doi/full/10.26434/chemrxiv-2025-0g4wn)

## License

MIT License — see `LICENSE`.

