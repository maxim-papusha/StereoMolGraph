# CondensedReactionGraph

### *class* stereomolgraph.CondensedReactionGraph(mol_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph) | None = None, frozen: bool = False)

Bases: [`MolGraph`](molgraph.md#stereomolgraph.MolGraph)

Graph representing a reaction. Atoms are nodes and (potentially changing)
bonds are edges. Every node has to have an attribute “atom_type” of type
Element. Edges can have an attribute “reaction” of type Change.
This is used to represent the change in connectivity during the reaction.

Two graphs are equal, iff. they are isomporhic and of the same type.

#### \_\_hash_\_() → int

Return hash(self).

* **Return type:**
  `int`

#### \_\_eq_\_(other: object) → bool

Return self==value.

* **Return type:**
  `bool`

#### add_bond(atom1: int, atom2: int, \*\*attr: Any)

Adds a bond between atom1 and atom2.

* **Parameters:**
  * **atom1** (`int`) – id of atom1
  * **atom2** (`int`) – id of atom2

#### set_bond_attribute(atom1: int, atom2: int, attr: str, value: Any)

sets the Attribute of the bond between Atom1 and Atom2.

* **Parameters:**
  * **atom1** (`int`) – Atom1
  * **atom2** (`int`) – Atom2
  * **attr** (`str`) – Attribute
  * **value** (`Any`) – Value

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

#### active_atoms(additional_layer: int = 0) → set[int]

Atoms involved in the reaction with additional layers of atoms
in the neighborhood.

* **Parameters:**
  **additional_layer** (`int`) – Number of additional layers of atoms to
  include, defaults to 0 (default: `0`)
* **Return type:**
  `set`[`int`]
* **Returns:**
  Atoms involved in the reaction

#### reactant(keep_attributes: bool = True) → [MolGraph](molgraph.md#stereomolgraph.MolGraph)

Reactant of the reaction

Creates the reactant of the reaction.
Formed bonds are not present in the reactant.

* **Parameters:**
  **keep_attributes** (`bool`) – attributes on atoms and bonds to be kept,
  defaults to True (default: `True`)
* **Return type:**
  [`MolGraph`](molgraph.md#stereomolgraph.MolGraph)
* **Returns:**
  Reactant of the reaction

#### product(keep_attributes: bool = True) → [MolGraph](molgraph.md#stereomolgraph.MolGraph)

Product of the reaction

Creates the product of the reaction.
Broken bonds are not present in the product.

* **Parameters:**
  **keep_attributes** (`bool`) – attributes on atoms and bonds to be kept,
  defaults to True (default: `True`)
* **Return type:**
  [`MolGraph`](molgraph.md#stereomolgraph.MolGraph)
* **Returns:**
  Product of the reaction

#### reverse_reaction() → Self

Creates the reaction in the opposite direction.

Broken bonds are turned into formed bonds and the other way around.

* **Return type:**
  `Self`
* **Returns:**
  Reversed reaction

#### *classmethod* from_graphs(reactant_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph), product_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph), ts_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph) | None = None) → Self

Creates a CondensedReactionGraph from reactant and product MolGraphs

CondensedReactionGraph  is constructed from bond changes from reactant
to the product. The atoms order and atom types of the reactant and
product have to be the same.

* **Parameters:**
  * **reactant_graph** ([`MolGraph`](molgraph.md#stereomolgraph.MolGraph)) – reactant of the reaction
  * **product_graph** ([`MolGraph`](molgraph.md#stereomolgraph.MolGraph)) – product of the reaction
* **Return type:**
  `Self`
* **Returns:**
  CondensedReactionGraph

#### *classmethod* from_geometries(reactant_geo: ~stereomolgraph.coords.GeometryProtocol, product_geo: ~stereomolgraph.coords.GeometryProtocol, ts_geo: ~stereomolgraph.coords.GeometryProtocol | None = None, switching_function: ~stereomolgraph.coords.BondsFromDistance = <stereomolgraph.coords.BondsFromDistance object>) → Self

Creates a CondensedReactionGraph from reactant
and product Geometries.

CondensedReactionGraph  is constructed from bond changes from reactant
to the product. The atoms order and atom types of the reactant and
product have to be the same. The switching function is used to
determine the connectivity of the atoms.

* **Parameters:**
  * **reactant_geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – geometry of the reactant
  * **product_geo** ([`GeometryProtocol`](../coords.md#stereomolgraph.coords.GeometryProtocol)) – geometry of the product
  * **switching_function** (`BondsFromDistance`) – function to define the connectivity
    from geometry,
    defaults to StepSwitchingFunction() (default: `<stereomolgraph.coords.BondsFromDistance object at 0x781c6c55f890>`)
* **Return type:**
  `Self`
* **Returns:**
  CondensedReactionGraph

#### \_\_init_\_(mol_graph: [MolGraph](molgraph.md#stereomolgraph.MolGraph) | None = None, frozen: bool = False)

#### *property* frozen *: bool*

Whether the graph is currently frozen (immutable).

#### freeze() → Self

Freeze the graph, making it immutable and hashable.

A frozen graph can be used in sets and as dict keys.
Mutation methods will raise `TypeError`.
Use [`copy()`](#stereomolgraph.CondensedReactionGraph.copy) to get a mutable copy.

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

Removes atom from graph.

* **Parameters:**
  **atom** (`int`) – Atom ID
* **Raises:**
  KeyError if atom is not in graph.

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

Returns a subgraph copy only containing the given atoms

* **Parameters:**
  **atoms** (`Iterable`[`int`]) – Iterable of atom ids to be
* **Return type:**
  `Self`
* **Returns:**
  Subgraph

#### copy(frozen: bool = False) → Self

* **Return type:**
  `Self`
* **Returns:**
  returns a copy of self

#### relabel_atoms(mapping: dict[int, int], copy: bool = True) → Self

Changes the atom labels according to mapping.

* **Parameters:**
  * **mapping** (`dict`[`int`, `int`]) – dict used for map old atom labels to new atom labels
  * **copy** (`bool`) – defines if the relabeling is done inplace or a new object
    should be created (default: `True`)
* **Return type:**
  `Self`
* **Returns:**
  this object (self) or a new instance of self._\_class_\_

#### *classmethod* from_rdmol(rdmol: Mol) → Self

Creates a StereoMolGraph from an RDKit Mol object.
Implicit Hydrogens are added to the graph.
Stereo information is conserved. Double bonds, aromatic bonds and
conjugated bonds are interpreted as planar. Atoms with 5 bonding
partners are assumed to be TrigonalBipyramidal and allow interchange
of the substituents (berry pseudorotation). Atoms with 6 bonding
partners are assumed to be octahedral and do not allow interchange of
the substituents.

* **Parameters:**
  **rdmol** (`Mol`) – RDKit Mol object
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

#### *classmethod* compose(mol_graphs: Iterable[[MolGraph](molgraph.md#stereomolgraph.MolGraph)]) → Self

Combines all graphs in the iterable into one. Duplicate nodes or edges
are overwritten, such that the resulting graph only contains one node
or edge with that name. Duplicate attributes of duplicate nodes or
edges are also overwritten in order of iteration.

* **Parameters:**
  **molgraphs** – Iterable of MolGraph that will be composed into a
  single MolGraph
* **Return type:**
  `Self`

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

Creates a graph of a molecule from a Geometry and a bond order matrix.

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
  Graph of Molecule

#### \_\_str_\_() → str

Return str(self).

* **Return type:**
  `str`

#### \_\_repr_\_()

Return repr(self).
