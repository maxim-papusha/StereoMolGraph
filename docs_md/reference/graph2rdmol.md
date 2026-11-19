# Graph2RDMol

### stereomolgraph.graph2rdmol.stereo_mol_graph_to_rdmol(graph: [StereoMolGraph](graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), generate_bond_orders=False, allow_charged_fragments=False, charge=0) → tuple[RWMol, dict[int, int]]

Creates a RDKit mol object using the connectivity of the mol graph.
Stereochemistry is added to the mol object.

* **Returns:**
  RDKit molecule
* **Return type:**
  Chem.rdchem.RWMol

### stereomolgraph.graph2rdmol.condensed_reaction_graph_to_rdmol(graph: [CondensedReactionGraph](graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph), generate_bond_orders: bool = False, allow_charged_fragments: bool = False, charge: int = 0) → tuple[RWMol, dict[int, int]]

Convert a CondensedReactionGraph to an RDKit molecule.

* **Parameters:**
  * **graph** ([`CondensedReactionGraph`](graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph)) – CondensedReactionGraph to convert
  * **generate_bond_orders** (`bool`) – If True, compute bond orders for the CRG (default: `False`)
  * **allow_charged_fragments** (`bool`) – If True, allow charged fragments (default: `False`)
  * **charge** (`int`) – Total charge of the molecule (default: `0`)
* **Return type:**
  `tuple`[`RWMol`, `dict`[`int`, `int`]]
* **Returns:**
  RDKit molecule and index-to-atom-id mapping

### stereomolgraph.graph2rdmol.stereo_condensed_reaction_graph_to_rdmol(graph: [StereoCondensedReactionGraph](graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph), generate_bond_orders: bool = False, allow_charged_fragments: bool = False, charge: int = 0) → tuple[RWMol, dict[int, int]]

Convert a StereoCondensedReactionGraph to an RDKit molecule.

The stereo changes are merged into a flat StereoMolGraph which is then
converted via [`stereo_mol_graph_to_rdmol()`](#stereomolgraph.graph2rdmol.stereo_mol_graph_to_rdmol).

* **Parameters:**
  * **graph** ([`StereoCondensedReactionGraph`](graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)) – StereoCondensedReactionGraph to convert
  * **generate_bond_orders** (`bool`) – If True, compute bond orders for the SCRG (default: `False`)
  * **allow_charged_fragments** (`bool`) – If True, allow charged fragments (default: `False`)
  * **charge** (`int`) – Total charge of the molecule (default: `0`)
* **Return type:**
  `tuple`[`RWMol`, `dict`[`int`, `int`]]
* **Returns:**
  RDKit molecule and index-to-atom-id mapping
