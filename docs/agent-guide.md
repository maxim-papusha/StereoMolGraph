# Using StereoMolGraph with coding agents

This is a task-oriented entry point for agent-generated Python code. The
[Sphinx API reference](reference/index.rst) provides detailed signatures.

## Install

\`\`\`bash
python -m pip install stereomolgraph
\`\`\`

## Main graph classes

Prefer importing from the package root:

\`\`\`python
from stereomolgraph import (
    MolGraph,
    StereoMolGraph,
    CondensedReactionGraph,
    StereoCondensedReactionGraph,
)
\`\`\`

- \`MolGraph\`: molecular connectivity without stereodescriptors.
- \`StereoMolGraph\`: molecular connectivity with local stereodescriptors.
- \`CondensedReactionGraph\`: graph representation of a reaction.
- \`StereoCondensedReactionGraph\`: reaction graphs including stereo changes.

## Convert an RDKit molecule

\`\`\`python
from rdkit import Chem
from stereomolgraph import StereoMolGraph

rdmol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](F)Cl"))
graph = StereoMolGraph.from_rdmol(rdmol)
print(graph.n_atoms)
print(graph.atom_stereo)
\`\`\`

Explicit hydrogens are required for bidirectional RDKit conversion.
\`StereoMolGraph.from_rdmol\` accepts \`stereo_complete\`, \`resonance\`, and
\`lone_pair_stereo\` keyword arguments. Consult its docstring for their meaning.

## Construct a connectivity graph

\`\`\`python
from stereomolgraph import MolGraph

graph = MolGraph()
graph.add_atom(0, "C")
graph.add_atom(1, "O")
graph.add_bond(0, 1)
\`\`\`

## Equality and mutability

- Graph equality tests graph isomorphism, not equality of atom identifiers.
  Stereo graphs additionally compare local stereochemistry.
- Freeze a graph before using it as a set element or dictionary key:

  \`\`\`python
  graph.freeze()
  unique_graphs = {graph}
  \`\`\`

## Scientific conventions and common pitfalls

- Bond orders, charge, and electronic state are not part of core graph identity.
  Graph equality is **not** a complete chemical identity comparison.
- Stereo descriptors represent relative local configurations. Parity values
  are **not** CIP \`R\`/\`S\` labels.
- Do not silently turn undefined or partially assigned stereochemistry into
  specified stereochemistry.
- Check public method signatures instead of assuming RDKit and StereoMolGraph
  APIs have the same semantics.

## Agent-readable documentation

The human-readable Sphinx documentation is the source of truth. It is
automatically converted into individual Markdown pages, \`llms.txt\` and
\`llms-full.txt\` using the \`sphinx-llm\` extension.
See the [repository agent instructions](https://github.com/maxim-papusha/StereoMolGraph/blob/doc_md/AGENTS.md).
