# StereoCondensedReactionGraph

### *class* stereomolgraph.StereoCondensedReactionGraph(mol_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph) | None = None)

Bases: [`StereoMolGraph`](stereomolgraph.md#stereomolgraph.StereoMolGraph), [`CondensedReactionGraph`](condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph)

`CondenedReactionGraph` with the ability to store stereochemistry
information for atoms and (potentially changing) bonds.

#### \_\_hash_\_() → int

Return hash(self).

* **Return type:**
  `int`

#### \_\_init_\_(mol_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph) | None = None)

#### \_\_eq_\_(other: object) → bool

Return self==value.

* **Return type:**
  `bool`

#### active_atoms(additional_layer: int = 0, stereo: bool = False) → set[int]

Atoms involved in the reaction with additional layers of atoms
in the neighborhood.

* **Parameters:**
  * **additional_layer** (`int`) – Number of additional layers of atoms to
    include, defaults to 0 (default: `0`)
  * **stereo** (`bool`) – Whether to include atoms involved in stereochemistry changes
    defaults to False (default: `False`)
* **Return type:**
  `set`[`int`]
* **Returns:**
  Atoms involved in the reaction

#### reactant(keep_attributes: bool = True) → [StereoMolGraph](stereomolgraph.md#stereomolgraph.StereoMolGraph)

Returns the reactant of the reaction

* **Parameters:**
  **keep_attributes** (`bool`) – If attributes should be kept , defaults to True (default: `True`)
* **Return type:**
  [`StereoMolGraph`](stereomolgraph.md#stereomolgraph.StereoMolGraph)
* **Returns:**
  reactant

#### product(keep_attributes: bool = True) → [StereoMolGraph](stereomolgraph.md#stereomolgraph.StereoMolGraph)

Returns the product of the reaction

* **Parameters:**
  **keep_attributes** (`bool`) – If attributes should be kept, defaults to True (default: `True`)
* **Return type:**
  [`StereoMolGraph`](stereomolgraph.md#stereomolgraph.StereoMolGraph)
* **Returns:**
  product

#### ts(infer_non_fleeting_stereo: bool = True) → [StereoMolGraph](stereomolgraph.md#stereomolgraph.StereoMolGraph)

Returns a StereoCondensedReactionGraph representing the transition state of
the reaction. Stereochemistry is taken from the fleeting stereo changes if
present, otherwise from the broken or formed stereo changes. If both broken and
formed stereo changes are present for a stereo element, the transition state
stereo is not defined and the stereo information is not included in the
transition state graph.

* **Return type:**
  [`StereoMolGraph`](stereomolgraph.md#stereomolgraph.StereoMolGraph)

#### reverse_reaction() → Self

Creates the reaction in the opposite direction.

Broken bonds and stereochemistry changes are turned into formed
and the other way around.

* **Return type:**
  `Self`
* **Returns:**
  Reversed reaction

#### copy(frozen: bool = False) → Self

* **Return type:**
  `Self`
* **Returns:**
  returns a copy of self

#### relabel_atoms(mapping: dict[int, int], copy: bool = True) → Self

Relabels the atoms of the graph and the chiral information accordingly

* **Parameters:**
  * **mapping** (`dict`[`int`, `int`]) – Mapping of old atom ids to new atom ids
  * **copy** (`bool`) – If the graph should be copied before relabeling,
    defaults to True (default: `True`)
* **Return type:**
  `Self`
* **Returns:**
  Returns the relabeled graph or None if copy is False

#### enantiomer() → Self

Creates the enantiomer of the StereoCondensedReactionGraph by inversion
of all chiral stereochemistries. The result can be identical to the
molecule itself if the molecule is not chiral.

* **Return type:**
  `Self`
* **Returns:**
  Enantiomer

#### *classmethod* compose(mol_graphs: Iterable[[MolGraph](molgraph.md#stereomolgraph.MolGraph)]) → Self

Creates a MolGraph object from a list of MolGraph objects

* **Parameters:**
  **mol_graphs** (`Iterable`[[`MolGraph`](molgraph.md#stereomolgraph.MolGraph)]) – list of MolGraph objects
* **Return type:**
  `Self`
* **Returns:**
  Returns Combined MolGraph

#### *property* frozen *: bool*

Whether the graph is currently frozen (immutable).

#### freeze() → Self

Freeze the graph, making it immutable and hashable.

A frozen graph can be used in sets and as dict keys.
Mutation methods will raise `TypeError`.
Use [`copy()`](#stereomolgraph.StereoCondensedReactionGraph.copy) to get a mutable copy.

* **Return type:**
  `Self`
* **Returns:**
  self (for chaining)

#### *property* atoms *: Collection[int]*

* **Returns:**
  Returns all atoms of the molecule

#### *property* atom_types *: tuple[int, ...]*

* **Returns:**
  Returns all atom types in the MolGraph

#### *property* atoms_with_attributes *: Mapping[int, dict[str, Any]]*

* **Returns:**
  Returns all atoms in the MolGraph with their attributes

#### *property* bonds *: Collection[frozenset[int]]*

* **Returns:**
  Returns all bonds in the MolGraph

#### *property* bonds_with_attributes *: Mapping[frozenset[int], dict[str, Any]]*

* **Returns:**
  Returns all bonds in the MolGraph with their attributes

#### *property* neighbors *: Mapping[int, frozenset[int]]*

* **Returns:**
  Returns all neighbors of the atoms in the MolGraph

#### *property* n_atoms *: int*

* **Returns:**
  Returns number of atoms in the MolGraph

#### has_atom(atom: int) → bool

Returns True if the molecules contains an atom with this id.

* **Parameters:**
  **atom** (`int`) – Atom
* **Return type:**
  `bool`
* **Returns:**
  value

#### add_atom(atom: int, atom_type: int | str, \*\*attr: Any)

Adds atom to the MolGraph

* **Parameters:**
  * **atom** (`int`) – Atom ID
  * **atom_type** (`int` | `str`) – Atom Type

#### remove_atom(atom: int)

Removes an atom from the graph and deletes all chiral information
associated with it

* **Parameters:**
  **atom** (`int`) – Atom

#### get_atom_attribute(atom: int, attr: str) → Any | None

Returns the value of the attribute of the atom or None if the atom does
not have this attribute.
Raises KeyError if atom is not in graph.

* **Parameters:**
  * **atom** (`int`) – Atom
  * **attr** (`str`) – Attribute
* **Raises:**
  **KeyError** – Atom not in graph
* **Return type:**
  `Optional`[`Any`]
* **Returns:**
  Returns the value of the attribute of the atom

#### get_atom_type(atom: int) → int

Returns the atom type of the atom.
Raises KeyError if atom is not in graph.

* **Parameters:**
  **atom** (`int`) – Atom
* **Raises:**
  **KeyError** – Atom not in graph
* **Return type:**
  `int`
* **Returns:**
  Returns the atom type of the atom

#### set_atom_attribute(atom: int, attr: str, value: Any)

sets the Value of the Attribute on Atom.
Raises KeyError if atom is not in graph.

* **Parameters:**
  * **atom** (`int`) – Atom
  * **attr** (`str`) – Attribute
  * **value** (`Any`) – Value
* **Raises:**
  * **KeyError** – Atom not in graph
  * **ValueError** – The attribute “atom_type” can only have values of
    type Element

#### delete_atom_attribute(atom: int, attr: str)

Deletes the Attribute of the Atom
Raises KeyError if attribute is not present.
Raises KeyError if atom is not in graph.

* **Parameters:**
  * **atom** (`int`) – Atom ID
  * **attr** (`str`) – Attribute
* **Raises:**
  **ValueError** – The attribute “atom_type” can not be deleted

#### get_atom_attributes(atom: int, attributes: Iterable[str] | None = None) → Mapping[str, Any]

Returns the attributes of the atom. If no attributes are given, all
attributes are returned.
Raises KeyError if atom is not in graph.

* **Parameters:**
  * **atom** (`int`) – Atom
  * **attributes** (`Optional`[`Iterable`[`str`]]) – Specific attributes to return (default: `None`)
* **Return type:**
  `Mapping`[`str`, `Any`]
* **Returns:**
  Returns all or just the chosen attributes of the atom

#### has_bond(atom1: int, atom2: int) → bool

Returns True if bond is in MolGraph.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2
* **Return type:**
  `bool`
* **Returns:**
  If the bond is in MolGraph

#### add_bond(atom1: int, atom2: int, \*\*attr: Any)

Adds a bond between atom1 and atom2.

* **Parameters:**
  * **atom1** (`int`) – id of atom1
  * **atom2** (`int`) – id of atom2

#### remove_bond(atom1: int, atom2: int)

Removes bond between Atom1 and Atom2.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2

#### get_bond_attribute(atom1: int, atom2: int, attr: str) → Any

Returns the value of the attribute of the bond between Atom1 and Atom2.
Raises KeyError if bond is not in graph.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2
  * **attr** (`str`) – Attribute
* **Return type:**
  `Any`
* **Returns:**
  Returns the value of the attribute of the bond
  between Atom1 and Atom2

#### set_bond_attribute(atom1: int, atom2: int, attr: str, value: Any)

sets the Attribute of the bond between Atom1 and Atom2.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2
  * **attr** (`str`) – Attribute
  * **value** (`Any`) – Value

#### delete_bond_attribute(atom1: int, atom2: int, attr: str)

Deletes the Attribute of the bond between Atom1 and Atom2

* **Parameters:**
  * **atom1** (`int`)
  * **atom2** (`int`) – Atom1
  * **attr** (`str`) – Attribute

#### get_bond_attributes(atom1: int, atom2: int, attributes: Iterable[str] | None = None) → Mapping[str, Any]

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2
  * **attributes** (`Optional`[`Iterable`[`str`]]) – Specific attributes to return (default: `None`)
* **Return type:**
  `Mapping`[`str`, `Any`]
* **Returns:**
  Returns chosen attributes of the bond between Atom1 and Atom2

#### get_atom_stereo(atom: int) → None | [AtomStereo](../stereodescriptors.md#stereomolgraph.stereodescriptors.AtomStereo)

Returns the stereo information of the atom if it exists else None.
Raises a ValueError if the atom is not in the graph.

* **Parameters:**
  * **atom** (`int`) – atom
  * **default** – Default value if no stereo information is found,
    defaults to None
* **Return type:**
  `None` | [`AtomStereo`](../stereodescriptors.md#stereomolgraph.stereodescriptors.AtomStereo)
* **Returns:**
  Stereo information of atom

#### set_atom_stereo(atom_stereo: [AtomStereo](../stereodescriptors.md#stereomolgraph.stereodescriptors.AtomStereo))

Adds stereo information to the graph

* **Parameters:**
  * **atom** – Atoms to be used for chiral information
  * **stereo** – Chiral information

#### delete_atom_stereo(atom: int)

Deletes stereo information from the graph

* **Parameters:**
  **atom** (`int`) – Atom to be used for stereo information

#### get_bond_stereo(bond: Iterable[int]) → None | [BondStereo](../stereodescriptors.md#stereomolgraph.stereodescriptors.BondStereo)

Gets the stereo information of the bond or None
if it does not exist.
Raises a ValueError if the bond s not in the graph.

* **Parameters:**
  **bond** (`Iterable`[`int`]) – Bond
* **Return type:**
  `None` | [`BondStereo`](../stereodescriptors.md#stereomolgraph.stereodescriptors.BondStereo)
* **Returns:**
  stereo information of bond

#### set_bond_stereo(bond_stereo: [BondStereo](../stereodescriptors.md#stereomolgraph.stereodescriptors.BondStereo))

Stets the stereo information of the bond

* **Parameters:**
  * **bond** – Bond
  * **bond_stereo** ([`BondStereo`](../stereodescriptors.md#stereomolgraph.stereodescriptors.BondStereo)) – Stereo information of the bond

#### delete_bond_stereo(bond: Iterable[int])

Deletes the stereo information of the bond

* **Parameters:**
  **bond** (`Iterable`[`int`]) – Bond

#### add_formed_bond(atom1: int, atom2: int, \*\*attr: Any)

Adds a bond between atom1 and atom2 with reaction attribute
set to FORMED.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2

#### add_broken_bond(atom1: int, atom2: int, \*\*attr: Any)

Adds a bond between atom1 and atom2 with reaction attribute
set to BROKEN.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2

#### add_fleeting_bond(atom1: int, atom2: int, \*\*attr: Any)

Adds a bond between atom1 and atom2 with reaction attribute
set to FLEETING.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2

#### get_formed_bonds() → set[frozenset[int]]

Returns all bonds that are formed during the reaction

* **Return type:**
  `set`[`frozenset`[`int`]]
* **Returns:**
  formed bonds

#### get_broken_bonds() → set[frozenset[int]]

Returns all bonds that are broken during the reaction

* **Return type:**
  `set`[`frozenset`[`int`]]
* **Returns:**
  broken bonds

#### get_fleeting_bonds() → set[frozenset[int]]

Returns all bonds that are fleeting during the reaction.

* **Return type:**
  `set`[`frozenset`[`int`]]
* **Returns:**
  fleeting bonds

#### bonded_to(atom: int) → frozenset[int]

Returns the atoms connected to the atom.

* **Parameters:**
  **atom** (`int`) – Id of the atom.
* **Return type:**
  `frozenset`[`int`]
* **Returns:**
  tuple of atoms connected to the atom.

#### connectivity_matrix() → ndarray[tuple[int, int], dtype[int8]]

Returns a connectivity matrix of the graph as a list of lists.
Order is the same as in self.atoms()
1 if nodes are connected, 0 if not.

* **Return type:**
  `ndarray`[`tuple`[`int`, `int`], `dtype`[`int8`]]
* **Returns:**
  Connectivity matrix as list of lists

#### node_connected_component(atom: int) → set[int]

* **Parameters:**
  **atom** (`int`) – atom id
* **Return type:**
  `set`[`int`]
* **Returns:**
  Returns the connected component that includes atom_id

#### connected_components() → list[set[int]]

* **Return type:**
  `list`[`set`[`int`]]
* **Returns:**
  Returns the connected components of the graph

#### subgraph(atoms: Iterable[int]) → Self

Returns a subgraph of the graph with the given atoms and the chiral
information accordingly. Stereo is retained when all referenced real
atoms remain; `None` lone-pair placeholders do not need to be selected.

* **Parameters:**
  **atoms** (`Iterable`[`int`]) – Atoms to be used for the subgraph
* **Return type:**
  `Self`
* **Returns:**
  Subgraph

#### *classmethod* from_rdmol(rdmol: Mol, , stereo_complete: bool = False, resonance: bool = True, lone_pair_stereo: bool = True) → Self

Creates a StereoMolGraph from an RDKit Mol object.
All hydrogens have to be explicit.
Stereo information is conserved for tetrahedral atoms and
double bonds.

* **Parameters:**
  * **rdmol** (`Mol`) – RDKit Mol object
  * **stereo_complete** (`bool`) – If True, we assume that the stereochemistry
    in the RDKit Mol is complete and all non chiral
    tetrahedral centers are set to an arbitrary
    configuration instead of None. (default: `False`)
  * **resonance** (`bool`) – If True, resonance structures are considered (default: `True`)
  * **lone_pair_stereo** (`bool`) – If True, lone pair stereo is considered (default: `True`)
* **Return type:**
  `Self`
* **Returns:**
  StereoMolGraph

#### bonds_from_bond_order_matrix(matrix: ndarray, threshold: float = 0.5, include_bond_order: bool = False)

Adds bonds the the graph based on bond orders from a matrix

* **Parameters:**
  * **matrix** (`ndarray`) – Bond order Matrix
  * **threshold** (`float`) – Threshold for bonds to be included as edges,
    defaults to 0.5 (default: `0.5`)
  * **include_bond_order** (`bool`) – If bond orders should be included as edge
    attributes, defaults to False (default: `False`)

#### *classmethod* from_atom_types_and_bond_order_matrix(atom_types: Sequence[int | str], matrix: ndarray, threshold: float = 0.5, include_bond_order: bool = False)

* **Parameters:**
  * **atom_types** (`Sequence`[`int` | `str`]) – list of atom types as integers or symbols,
    must correspond to the matrix
  * **matrix** (`ndarray`) – np.matrix of bond orders or connectivities ([0..1])
  * **threshold** (`float`) – Threshold for bonds to be included as edges,
    defaults to 0.5 (default: `0.5`)
  * **include_bond_order** (`bool`) – If bond orders should be included as edge
    attributes, defaults to False (default: `False`)
* **Returns:**
  Returns MolGraph

#### *classmethod* from_geometry_and_bond_order_matrix(geo: [GeometryProtocol](../coords.md#stereomolgraph.coords.GeometryProtocol), matrix: ndarray, threshold: float = 0.5, include_bond_order: bool = False) → Self

Creates a CiralMolGraph object from a Geometry and a bond order matrix

* **Parameters:**
  * **geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – Geometry
  * **matrix** (`ndarray`) – Bond order matrix
  * **threshold** (`float`) – Threshold for bonds to be included as edges,
    defaults to 0.5 (default: `0.5`)
  * **include_bond_order** (`bool`) – If bond orders should be included as edge
    attributes, defaults to False (default: `False`)
* **Return type:**
  `Self`
* **Returns:**
  Returns MolGraph

#### *classmethod* from_geometry(geo: ~stereomolgraph.coords.GeometryProtocol, switching_function: ~stereomolgraph.coords.BondsFromDistance = <stereomolgraph.coords.BondsFromDistance object>) → Self

Creates a StereoMolGraph object from a Geometry and a switching
function. Uses the Default switching function if none are given.

* **Parameters:**
  * **geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – Geometry
  * **switching_function** (`BondsFromDistance`) – Function to determine if two atoms are
    connected (default: `<stereomolgraph.coords.BondsFromDistance object at 0x79a9b7d06190>`)
* **Return type:**
  `Self`
* **Returns:**
  StereoMolGraph of molecule

#### *classmethod* from_graphs(reactant_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph), product_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph), ts_graph: None | [MolGraph](molgraph.md#stereomolgraph.MolGraph) = None) → Self

Creates a StereoCondensedReactionGraph from reactant and product
StereoMolGraphs.

StereoCondensedReactionGraph  is constructed from bond changes from
reactant to the product. The atoms order and atom types of the reactant
and product have to be the same.

* **Parameters:**
  * **reactant_graph** ([`MolGraph`](molgraph.md#stereomolgraph.MolGraph)) – reactant of the reaction
  * **product_graph** ([`MolGraph`](molgraph.md#stereomolgraph.MolGraph)) – product of the reaction
* **Return type:**
  `Self`
* **Returns:**
  StereoCondensedReactionGraph

#### \_\_str_\_() → str

Return str(self).

* **Return type:**
  `str`

#### \_\_repr_\_()

Return repr(self).

#### *classmethod* from_geometries(reactant_geo: ~stereomolgraph.coords.GeometryProtocol, product_geo: ~stereomolgraph.coords.GeometryProtocol, ts_geo: None | ~stereomolgraph.coords.GeometryProtocol = None, switching_function: ~stereomolgraph.coords.BondsFromDistance = <stereomolgraph.coords.BondsFromDistance object>) → Self

Creates a StereoCondensedReactionGraph from reactant, product and
transition state Geometries.

StereoCondensedReactionGraph  is constructed from bond changes from
reactant to the product. The atoms order and atom types of the reactant
and product have to be the same. The switching function is used to
determine the connectivity of the atoms. Only the stereo information
is taken from the transition state geometry.

* **Parameters:**
  * **reactant_geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – geometry of the reactant
  * **product_geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – geometry of the product
  * **ts_geo** (`None` | [`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – geometry of the transition state (default: `None`)
  * **switching_function** (`BondsFromDistance`) – function to define the connectivity from
    geometry,
    defaults to StepSwitchingFunction() (default: `<stereomolgraph.coords.BondsFromDistance object at 0x79a9b7d062c0>`)
* **Return type:**
  `Self`
* **Returns:**
  CondensedReactionGraph
