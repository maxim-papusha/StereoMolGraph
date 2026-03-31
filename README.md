![Logo](https://raw.githubusercontent.com/maxim-papusha/StereoMolGraph/main/docs/_static/img/logo_smg.png)

# StereoMolGraph #

[![PyPI](https://img.shields.io/pypi/v/StereoMolGraph?style=flat-square&logo=pypi&logoColor=white&color=3775A9)](https://pypi.org/project/StereoMolGraph/)
[![Python](https://img.shields.io/pypi/pyversions/StereoMolGraph?style=flat-square&logo=python&logoColor=white&color=3776AB&label=Python)](https://pypi.org/project/StereoMolGraph/)
[![License: MIT](https://img.shields.io/badge/License-MIT-d4a017?style=flat-square&logo=opensourceinitiative&logoColor=white)](https://opensource.org/licenses/MIT)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/maxim-papusha/StereoMolGraph/run_unit_test.yaml?branch=main&style=flat-square&label=tests)](https://github.com/maxim-papusha/StereoMolGraph/actions/workflows/run_unit_test.yaml)

[![Documentation](https://img.shields.io/badge/Documentation-docs-4C6A92?style=flat-square&logo=readthedocs&logoColor=white)](https://stereomolgraph.readthedocs.io)
[![GitHub](https://img.shields.io/badge/GitHub-repository-2F3E46?style=flat-square&logo=github&logoColor=white)](https://github.com/maxim-papusha/StereoMolGraph)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5c02523-ffcc00?style=flat-square)](https://doi.org/10.1021/acs.jcim.5c02523)

StereoMolGraph (SMG) is a library for graph representation of molecules and reactions with a focus on Stereochemistry. It provides:

- Graph types for molecules and reactions (with/without stereo and stereo changes)
- Includes non tetrahedral stereocenters and changing stereochemistry in reactions
- Fast approximate hashing via Weisfeiler–Lehman color refinement
- Robust equality/isomorphism via a VF2++-style algorithm extended for stereochemistry and reactions
- Bidirectional conversion from / to RDKit
- Construction from 3D coordinates with automatic local stereo inference


## Design philosophy

- Unopinionated about bond orders, charge and electronic state
- SMG focuses just on the connectivity and stereochemistry. 
- Stereochemistry describes relative spatial arrangement. No absolute stereochemistry.
- Transparent: Simple 2D visualization in IPython notebooks

## RDKit interoperability notes

- Hydrogens must be explicit for bidirectional conversion.
- Supports tetrahedral and non tetrahedral stereochemistry during conversion.
- Bond orders, charges, unpaired electrons and other properties are not used!

## Installation

Install from PyPI:

```bash
pip install stereomolgraph
```

## Feedback and support

Bug reports are welcome — please open an issue on GitHub.
Technical questions and feature requests can be added as discussions.

Issues:
- Bugreports: https://github.com/maxim-papusha/StereoMolGraph/issues

Discussions:
- Q&A: https://github.com/maxim-papusha/StereoMolGraph/discussions/new?category=q-a
- Ideas: https://github.com/maxim-papusha/StereoMolGraph/discussions/new?category=ideas


## Citation

If you use StereoMolGraph in your work, please cite the corresponding paper:

[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5c02523-ffcc00?style=flat-square)](https://doi.org/10.1021/acs.jcim.5c02523)

## License

MIT License — see `LICENSE`.

