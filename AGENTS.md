# StereoMolGraph — instructions for coding agents

## Purpose
StereoMolGraph provides graph representations of molecules and reactions,
including local stereochemical descriptors, stereochemistry-aware isomorphism,
and symmetry analysis.

## Read first
- [Agent usage guide](docs_md/agent-guide.md)
- [Generated Markdown reference](docs_md/llms.txt)
- [Browsable API documentation](https://stereomolgraph.readthedocs.io/)

## Public Python API
Prefer importing from the package root:
\`\`\`python
from stereomolgraph import (
    MolGraph,
    StereoMolGraph,
    CondensedReactionGraph,
    StereoCondensedReactionGraph,
)
\`\`\`

## Chemical invariants
- Core molecular graphs describe connectivity and local stereochemistry.
- Bond orders and formal charges are not part of core graph identity.
- Stereo parity is a local relative descriptor, not a CIP R/S assignment.
- Do not assign undefined stereochemistry implicitly.
- RDKit conversion requires explicit hydrogens for faithful round trips.
- Graph equality is based on isomorphism; atom IDs do not establish identity.
- Mutable graphs are unhashable. Call .freeze() before set/dict-key use.

## Contribution rules
- Preserve existing graph identity and stereochemical equivalence semantics.
- Prefer small public APIs and task-based executable examples.
- Read docstrings and tests before changing algorithms.
- Test with pytest and lint with ruff check .
- The canonical documentation is in docs/ and source docstrings.
- docs_md/ and the root llms.txt are generated; do not edit them manually.
- GitHub Actions rebuilds Markdown on pushes to the doc_md branch.

## Local documentation build
\`\`\`bash
python -m pip install -e '.[docs]'
python -m sphinx -b html docs docs/_build/html
python scripts/export_markdown_docs.py
\`\`\`
