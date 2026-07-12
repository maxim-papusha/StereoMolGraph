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
- Fast hashing using Circular Stereo Hash
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

Issues:
- Bugreports: https://github.com/maxim-papusha/StereoMolGraph/issues


## Citation

If you use **StereoMolGraph** in your research, please cite the main library paper:

> **M. Papusha and K. Leonhard**,
> “StereoMolGraph: Stereochemistry-Aware Molecular and Reaction Graphs,”
> *Journal of Chemical Information and Modeling* **2026**, *66* (7), 3830–3839.

[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5c02523-ffcc00?style=flat-square)](https://doi.org/10.1021/acs.jcim.5c02523)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-StereoMolGraph)

## Publications using StereoMolGraph

The following publications use StereoMolGraph. The corresponding repositories contain experimental code, notebooks, data, and reproducibility material.

### Circular Stereo Algorithm and Fingerprint

> **M. Papusha and K. Leonhard**,
> “Circular Stereo Algorithm and Fingerprint for Chiral Resonance Invariant Molecular Representation,”
> *ChemRxiv* **2026**, preprint.

[![Preprint](https://img.shields.io/badge/Preprint-ChemRxiv-B31B1B?style=flat-square)](https://doi.org/10.26434/chemrxiv.15002723/v1)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-CircularStereoAlgorithm)

### Symmetry Numbers

> **M. Papusha, A. V. Copan, B. Rotavera, and K. Leonhard**,
> “Symmetry Numbers: A Flexible Approach for Molecules and Transition States,”
> manuscript submitted for publication, **2026**.

[![Status](https://img.shields.io/badge/Status-Submitted-6B7280?style=flat-square)](#)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-SymmetryNumbers)

## BibTeX

<details>
<summary><strong>Show all BibTeX references</strong></summary>

```bibtex
@article{Papusha2026StereoMolGraph,
  author  = {Papusha, Maxim and Leonhard, Kai},
  title   = {{StereoMolGraph}: Stereochemistry-Aware Molecular and Reaction Graphs},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2026},
  volume  = {66},
  number  = {7},
  pages   = {3830--3839},
  doi     = {10.1021/acs.jcim.5c02523},
  url     = {https://doi.org/10.1021/acs.jcim.5c02523}
}

@article{Papusha2026CircularStereoAlgorithm,
  author  = {Papusha, Maxim and Leonhard, Kai},
  title   = {Circular Stereo Algorithm and Fingerprint for Chiral Resonance Invariant Molecular Representation},
  journal = {ChemRxiv},
  year    = {2026},
  note    = {Preprint},
  doi     = {10.26434/chemrxiv.15002723/v1},
  url     = {https://doi.org/10.26434/chemrxiv.15002723/v1}
}

@unpublished{Papusha2026SymmetryNumbers,
  author = {Papusha, Maxim and Copan, Andreas V. and Rotavera, Brandon and Leonhard, Kai},
  title  = {Symmetry Numbers: A Flexible Approach for Molecules and Transition States},
  year   = {2026},
  note   = {Manuscript submitted for publication},
  url    = {https://github.com/maxim-papusha/Experiments-SymmetryNumbers}
}
```

</details>



## License

MIT License — see `LICENSE`.

