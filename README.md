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

Bug reports, feature requests, and questions are welcome through
[GitHub Issues](https://github.com/maxim-papusha/StereoMolGraph/issues).


## Citations

If you use **StereoMolGraph**, please cite the relevant publication(s):

1. M. Papusha and K. Leonhard, “**StereoMolGraph**: Stereochemistry-Aware Molecular and Reaction Graphs,” *J. Chem. Inf. Model.*, 2026, *66*, 3830–3839.
   [![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5c02523-ffcc00?style=flat-square)](https://doi.org/10.1021/acs.jcim.5c02523)
   [![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-StereoMolGraph)

2. M. Papusha and K. Leonhard, “**Circular Stereo Algorithm** and Fingerprint for Chiral Resonance Invariant Molecular Representation,” *ChemRxiv*, 2026, preprint.
   [![DOI](https://img.shields.io/badge/DOI-10.26434%2Fchemrxiv.15002723%2Fv1-ffcc00?style=flat-square)](https://doi.org/10.26434/chemrxiv.15002723/v1)
   [![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-CircularStereoAlgorithm)

3. M. Papusha, A. V. Copan, B. Rotavera, and K. Leonhard, “**Symmetry Numbers**: A Flexible Approach for Molecules and Transition States,” submitted manuscript, 2026.
   ![DOI](https://img.shields.io/badge/DOI-pending-lightgrey?style=flat-square)
   [![Reproducibility](https://img.shields.io/badge/Reproducibility-Code-181717?style=flat-square\&logo=github\&logoColor=white)](https://github.com/maxim-papusha/Experiments-SymmetryNumbers)

<details>
<summary><strong>BibTeX entries</strong></summary>

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
}

@article{Papusha2026CircularStereoAlgorithm,
  author  = {Papusha, Maxim and Leonhard, Kai},
  title   = {Circular Stereo Algorithm and Fingerprint for Chiral Resonance
             Invariant Molecular Representation},
  journal = {ChemRxiv},
  year    = {2026},
  note    = {Preprint},
  doi     = {10.26434/chemrxiv.15002723/v1},
}

@unpublished{Papusha2026SymmetryNumbers,
  author = {Papusha, Maxim and Copan, Andreas V. and
            Rotavera, Brandon and Leonhard, Kai},
  title  = {Symmetry Numbers: A Flexible Approach for Molecules
            and Transition States},
  year   = {2026},
  note   = {Manuscript submitted for publication},
}
```

</details>


## License

MIT License — see `LICENSE`.

