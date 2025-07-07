"""Graph theory

This module contains classes for representing graphs and performing graph
theory operations on them. It is meant as general, application-agnostic code
and should not contain graph classes designed for a single purpose only.
"""

from __future__ import annotations

import warnings
from collections import Counter, defaultdict, deque
from copy import deepcopy
from enum import Enum
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np

from stereomolgraph import PERIODIC_TABLE, Element
from stereomolgraph.cartesian import are_planar, BondsFromDistance
from stereomolgraph.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.color_refine import color_refine_mg
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    BondStereo,
    Tetrahedral,
    TrigonalBipyramidal,
    Octahedral,
    PlanarBond,
    Stereo,
    SquarePlanar)


from stereomolgraph.graph2rdmol import (
    _mol_graph_to_rdmol,
    _stereo_mol_graph_to_rdmol,
    _set_crg_bond_orders
)
from stereomolgraph.rdmol2graph import (
    mol_graph_from_rdmol, 
    stereo_mol_graph_from_rdmol,
    )

if TYPE_CHECKING:

    import sys
    from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
    from typing import Any, Optional, TypeAlias
    
    from rdkit import Chem
    import scipy.sparse  # type: ignore

    from stereomolgraph.cartesian import Geometry
       
    # Self is included in typing from 3.11
    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self


AtomId: TypeAlias = int
    
class Bond(tuple[int, int]):

    def __new__(cls, atoms: Iterable[int]):
        ret = super().__new__(cls, sorted(atoms))
        if len(ret) != 2 or ret[0] == ret[1]:
            raise ValueError("A bond has to connect two distinct atoms")
        return ret


class MolGraph:
    """
    Graph representing a molecular entity. Nodes represent atoms and edges
    represent bonds. All nodes have an `atom_type` attribute of type `Element`.
    The node ids should be integers. The graph is considered equal to another
    graph, iff. they are isomorphic and of the same type.
    """
    __slots__ = ("_atom_attrs", "_neighbors", "_bond_attrs")

    _atom_attrs: dict[AtomId, dict[str, Any]]
    _neighbors: dict[AtomId, set[AtomId]]
    _bond_attrs: dict[Bond, dict[str, Any]]

    def __init__(self, mol_graph: Optional[MolGraph] = None):
        if mol_graph is not None:
            if not isinstance(mol_graph, MolGraph):
                raise ValueError("Input has to be of type MolGraph")
            self._atom_attrs = deepcopy(mol_graph._atom_attrs)
            self._neighbors = deepcopy(mol_graph._neighbors)
            self._bond_attrs = deepcopy(mol_graph._bond_attrs)
        else:
            self._atom_attrs = defaultdict(dict)
            self._neighbors = defaultdict(set)
            self._bond_attrs = defaultdict(dict)

    @property
    def atoms(
        self,
    ) -> Iterable[AtomId]:
        """
        :return: Returns all atoms of the molecule
        """
        return self._atom_attrs.keys()

    @property
    def atom_types(
        self,
    ) -> tuple[Element, ...]:
        """
        :return: Returns all atom types in the MolGraph
        """
        return tuple([v["atom_type"] for v in self._atom_attrs.values()])

    @property
    def atoms_with_attributes(self) -> Mapping[AtomId, dict[str, Any]]:
        """
        :return: Returns all atoms in the MolGraph with their attributes
        """
        return MappingProxyType(self._atom_attrs)

    @property
    def bonds(
        self,
    ) -> Iterable[Bond]:
        """
        :return: Returns all bonds in the MolGraph
        """
        return self._bond_attrs.keys()

    @property
    def bonds_with_attributes(
        self,
    ) -> Mapping[Bond, dict[str, Any]]:
        """
        :return: Returns all bonds in the MolGraph with their attributes
        """
        return MappingProxyType(self._bond_attrs)

    @property
    def n_atoms(
        self,
    ) -> int:
        """
        :return: Returns number of atoms in the MolGraph
        """
        return len(self._atom_attrs)

    def __len__(self) -> int:
        return len(self._atom_attrs)

    def has_atom(self, atom: int) -> bool:
        """Returns True if the molecules contains an atom with this id.

        :param atom: Atom
        :return: value
        """
        return atom in self._atom_attrs

    def add_atom(
        self, atom: AtomId, atom_type: int | str | Element, **attr
    ):
        """Adds atom to the MolGraph

        :param atom: Atom ID
        :param atom_type: Atom Type
        """
        atom_type = PERIODIC_TABLE[atom_type]

        self._atom_attrs[atom] = {"atom_type": atom_type, **attr}

    def remove_atom(self, atom: AtomId):
        """Removes atom from graph.
        
        :param atom: Atom ID
        :raises: KeyError if atom is not in graph.
        """
        del self._atom_attrs[atom]
        if nbr := self._neighbors.pop(atom, None):
            for n in nbr:
                self._neighbors[n].remove(atom)

    def get_atom_attribute(
        self, atom: AtomId, attr: str
    ) -> Optional[Any]:
        """
        Returns the value of the attribute of the atom or None if the atom does
        not have this attribute.
        Raises KeyError if atom is not in graph.

        :param atom: Atom
        :param attr: Attribute
        :raises KeyError: Atom not in graph
        :return: Returns the value of the attribute of the atom
        """
        return self._atom_attrs[atom].get(attr, None)

    def get_atom_type(self, atom: AtomId) -> Element:
        """
        Returns the atom type of the atom.
        Raises KeyError if atom is not in graph.

        :param atom: Atom
        :raises KeyError: Atom not in graph
        :return: Returns the atom type of the atom
        """
        return self._atom_attrs[atom]["atom_type"]

    def set_atom_attribute(self, atom: AtomId, attr: str, value: Any):
        """
        sets the Value of the Attribute on Atom.
        Raises KeyError if atom is not in graph.

        :param atom: Atom
        :param attr: Attribute
        :param value: Value
        :raises KeyError: Atom not in graph
        :raises ValueError: The attribute "atom_type" can only have values of
                            type Element
        """
        if attr == "atom_type":
            try:
                value = PERIODIC_TABLE[value]
            except KeyError:
                raise ValueError(
                    f"'{value}' can not be used as atom_type for "
                    f"{self.__class__.__name__}"
                )
        self._atom_attrs[atom][attr] = value

    def delete_atom_attribute(self, atom: AtomId, attr: str):
        """
        Deletes the Attribute of the Atom
        Raises KeyError if attribute is not present.
        Raises KeyError if atom is not in graph.

        :param atom: Atom ID
        :param attr: Attribute
        :raises ValueError: The attribute "atom_type" can not be deleted
        """
        if attr == "atom_type":
            raise ValueError("atom_type can not be deleted")
        else:
            self._atom_attrs[atom].pop(attr)

    def get_atom_attributes(
        self, atom: AtomId, attributes: Optional[Iterable[str]] = None
    ) -> Mapping[str, Any]:
        """
        Returns the attributes of the atom. If no attributes are given, all
        attributes are returned.
        Raises KeyError if atom is not in graph.

        :param atom: Atom
        :param attributes: Specific attributes to return
        :return: Returns all or just the chosen attributes of the atom
        """
        if attributes is None:
            return MappingProxyType(self._atom_attrs[atom])
        else:
            return {attr: self._atom_attrs[atom][attr] for attr in attributes}

    def has_bond(self, atom1: AtomId, atom2: AtomId) -> bool:
        """Returns True if bond is in MolGraph.

        :param atom1: Atom1
        :param atom2: Atom2
        :return: If the bond is in MolGraph
        """
        return Bond((atom1, atom2)) in self._bond_attrs

    def add_bond(self, atom1: AtomId, atom2: AtomId, **attr):
        """Adds bond between Atom1 and Atom2.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        if atom1 not in self.atoms or atom2 not in self.atoms:
            raise ValueError("Atoms not in Graph")
        bond = Bond({atom1, atom2})
        self._neighbors[atom1].add(atom2)
        self._neighbors[atom2].add(atom1)
        self._bond_attrs[bond] = attr

    def remove_bond(self, atom1: AtomId, atom2: AtomId):
        """
        Removes bond between Atom1 and Atom2.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        bond = Bond((atom1, atom2))
        del self._bond_attrs[bond]
        self._neighbors[atom1].remove(atom2)
        self._neighbors[atom2].remove(atom1)

    def get_bond_attribute(
        self, atom1: AtomId, atom2: AtomId, attr: str,
    ) -> Any:
        """
        Returns the value of the attribute of the bond between Atom1 and Atom2.
        Raises KeyError if bond is not in graph.

        :param atom1: Atom1
        :param atom2: Atom2
        :param attr: Attribute
        :return: Returns the value of the attribute of the bond
                 between Atom1 and Atom2
        """
        bond = Bond((atom1, atom2))
        if bond in self._bond_attrs:
            return self._bond_attrs[bond].get(attr, None)
        else:
            raise ValueError(f"No Bond between {atom1} and {atom2}")

    def set_bond_attribute(
        self, atom1: AtomId, atom2: AtomId, attr: str, value: Any
    ):
        """
        sets the Attribute of the bond between Atom1 and Atom2.
        The Attribute "bond_order" can only have numerical values.
        Raises KeyError if bond is not in graph.

        :param atom1: Atom1
        :param atom2: Atom2
        :param attr: Attribute
        :param value: Value
        """
        bond = Bond((atom1, atom2))
        if bond in self._bond_attrs:
            self._bond_attrs[bond][attr] = value
        else:
            raise ValueError(f"No Bond between {atom1} and {atom2}")

    def delete_bond_attribute(self, atom1: AtomId, atom2: AtomId,
                              attr: str):
        """
        Deletes the Attribute of the bond between Atom1 and Atom2
        
        :param atom1:
        :param atom2: Atom1
        :param attr: Attribute
        """
        self._bond_attrs[Bond((atom1, atom2))].pop(attr)

    def get_bond_attributes(
        self,
        atom1: AtomId,
        atom2: AtomId,
        attributes: Optional[Iterable[str]] = None,
    ) -> Mapping[str, Any]:
        """
        :param atom1: Atom1
        :param atom2: Atom2
        :param attributes: Specific attributes to return
        :return: Returns chosen attributes of the bond between Atom1 and Atom2
        """
        bond = Bond((atom1, atom2))
        if attributes is None:
            return MappingProxyType(self._bond_attrs[bond])
        else:
            return {attr: val for attr, val in self._bond_attrs[bond].items()}

    def bonded_to(self, atom: int) -> frozenset[int]:
        """
        Returns the atoms connected to the atom.

        :param atom: Id of the atom.
        :return: tuple of atoms connected to the atom.
        """
        return frozenset(self._neighbors[atom])

    def connectivity_matrix(self) -> list[list[int]]:
        """
        Returns a connectivity matrix of the graph as a list of lists. Order is the same as
        in self.atoms()
        1 if nodes are connected, 0 if not.

        :return: Connectivity matrix as list of lists
        """
        matrix = [[0] * self.n_atoms for _ in range(self.n_atoms)]
        for atom1, atom2 in self.bonds:
            matrix[atom1][atom2] = 1
            matrix[atom2][atom1] = 1
        return matrix

    def _to_rdmol(
        self, generate_bond_orders=False, charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        return _mol_graph_to_rdmol(self,
                                   generate_bond_orders=generate_bond_orders,
                                   charge=charge)

    def to_rdmol(
        self, generate_bond_orders=False, charge=0
    ) -> Chem.rdchem.Mol:

        mol, _ = self._to_rdmol(
            generate_bond_orders=generate_bond_orders, charge=charge
        )
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0, strict=True)
        return mol

    @classmethod
    def from_rdmol(cls, rdmol, use_atom_map_number=False) -> Self:
        """
        Creates a StereoMolGraph from an RDKit Mol object.
        Implicit Hydrogens are added to the graph.
        Stereo information is conserved. Double bonds, aromatic bonds and
        conjugated bonds are interpreted as planar. Atoms with 5 bonding
        partners are assumed to be TrigonalBipyramidal and allow interchange
        of the substituents (berry pseudorotation). Atoms with 6 bonding
        partners are assumed to be octahedral and do not allow interchange of
        the substituents.

        :param rdmol: RDKit Mol object
        :param use_atom_map_number: If the atom map number should be used
                                    instead of the atom index
        :return: StereoMolGraph
        """
        return mol_graph_from_rdmol(cls, rdmol,
                                    use_atom_map_number=use_atom_map_number)

    def relabel_atoms(
        self, mapping: dict[int, int], copy: bool = True
    ) -> Self:
        """Changes the atom labels according to mapping.

        :param mapping: dict used for map old atom labels to new atom labels
        :param copy: defines if the relabeling is done inplace or a new object
                     should be created
        :return: this object (self) or a new instance of self.__class__
        """
        atom_attrs = {
            mapping.get(atom, atom): attrs
            for atom, attrs in self._atom_attrs.items()
        }
        neighbors = {
            mapping.get(atom, atom): {mapping.get(n, n) for n in neighbors}
            for atom, neighbors in self._neighbors.items()
        }

        bond_attrs = {
            Bond({mapping.get(atom, atom) for atom in bond}): attrs
            for bond, attrs in self._bond_attrs.items()
        }
        if copy is True:
            new_graph = self.__class__()
        elif copy is False:
            new_graph = self

        new_graph._atom_attrs = atom_attrs
        new_graph._neighbors = neighbors
        new_graph._bond_attrs = bond_attrs
        return new_graph

    def node_connected_component(self, atom: int) -> set[int]:
        """
        :param atom: atom id
        :return: Returns the connected component that includes atom_id
        """
        visited = set()
        for layer in self.bfs_layers(atom):
            for node in layer:
                if node not in visited:
                    visited.add(node)
        return visited

    def connected_components(self) -> list[set[int]]:
        """
        :return: Returns the connected components of the graph
        """
        visited = set()
        components = []

        for atom in self.atoms:
            if atom not in visited:
                component = self.node_connected_component(atom)
                components.append(component)
                visited.update(component)

        return components

    def subgraph(self, atoms: Iterable[AtomId]) -> Self:
        """
        Returns a subgraph copy only containing the given atoms

        :param atoms: Iterable of atom ids to be
        :return: Subgraph
        """
        new_atoms = set(atoms)
        atom_attrs = {atom: self._atom_attrs[atom] for atom in atoms}
        bond_attrs = {
            bond: attrs
            for bond, attrs in self._bond_attrs.items()
            if new_atoms.issuperset(bond)
        }
        neighbors = {
            atom: {n for n in self._neighbors[atom] if n in new_atoms}
            for atom in new_atoms
        }
        new_graph = self.__class__()
        new_graph._atom_attrs = atom_attrs
        new_graph._neighbors = neighbors
        new_graph._bond_attrs = bond_attrs
        return new_graph

    def copy(self) -> Self:
        """
        :return: returns a copy of self
        """
        return deepcopy(self)

    def bonds_from_bond_order_matrix(
        self,
        matrix: np.ndarray | scipy.sparse.sparray,
        threshold: float = 0.5,
        include_bond_order: bool = False,
    ):
        """
        Adds bonds the the graph based on bond orders from a matrix

        :param matrix: Bond order Matrix
        :param threshold: Threshold for bonds to be included as edges, 
                          defaults to 0.5
        :param include_bond_order: If bond orders should be included as edge
                                   attributes, defaults to False
        """

        if not np.shape(matrix) == (len(self), len(self)):
            raise ValueError(
                "Matrix has the wrong shape. shape of matrix is "
                f"{np.shape(matrix)}, but {len(self), len(self)} "
                "expected"
            )

        bonds = (matrix > threshold).nonzero()

        for i, j in zip(*bonds):
            if include_bond_order:
                self.add_bond(int(i), int(j), bond_order=matrix[i, j])
            else:
                self.add_bond(int(i), int(j))

    @classmethod
    def from_composed_molgraphs(cls, mol_graphs: Iterable[Self]) -> Self:
        """
        Combines all graphs in the iterable into one. Duplicate nodes or edges
        are overwritten, such that the resulting graph only contains one node
        or edge with that name. Duplicate attributes of duplicate nodes or
        edges are also overwritten in order of iteration.

        :param molgraphs: Iterable of MolGraph that will be composed into a
            single MolGraph
        """
        new_graph = cls()
        for mol_graph in mol_graphs:
            new_graph._atom_attrs.update(mol_graph._atom_attrs)
            new_graph._bond_attrs.update(mol_graph._bond_attrs)

            for atom, neighbors in mol_graph._neighbors.items():
                new_graph._neighbors[atom].update(neighbors)

        return new_graph

    @classmethod
    def from_atom_types_and_bond_order_matrix(
        cls,
        atom_types: Sequence[int | Element | str],
        matrix: np.ndarray,
        threshold=0.5,
        include_bond_order=False,
    ):
        """

        :param atom_types: list of atom types as integers or symbols,
                           must correspond to the matrix
        :param matrix: np.matrix of bond orders or connectivities ([0..1])
        :param threshold: Threshold for bonds to be included as edges,
                          defaults to 0.5
        :param include_bond_order: If bond orders should be included as edge
                                   attributes, defaults to False
        :return: Returns MolGraph
        """
        if not len(atom_types) == np.shape(matrix)[0] == np.shape(matrix)[1]:
            raise ValueError(
                "atom_types and matrix have to have the same length"
            )
        new_mol_graph = cls()

        for i, atom_type in enumerate(atom_types):
            new_mol_graph.add_atom(i, atom_type=atom_type)

        x_ids, y_ids = np.triu_indices(matrix.shape[0], k=1)

        # Iterate over the upper triangular matrix excluding the diagonal
        for x_id, y_id in zip(x_ids, y_ids):
            if (value := matrix[x_id, y_id]) >= threshold:
                if include_bond_order is False:
                    new_mol_graph.add_bond(int(x_id), int(y_id))
                elif include_bond_order is True:
                    new_mol_graph.add_bond(
                        int(x_id), int(y_id), bond_order=value
                    )

        return new_mol_graph

    @classmethod
    def from_atom_types_and_bond_order_sparse_array(
        cls,
        atom_types: Sequence[int | Element | str],
        sp_arr: scipy.sparse.sparray,
        threshold: float = 0.5,
        include_bond_order=False,
    ):
        """

        :param atom_types: list of atom types as integers or symbols,
                           must correspond to the array
        :param matrix: scipy.sparse.sparray of bond orders or connectivities
                       ([0..1])
        :param threshold: Threshold for bonds to be included as edges,
                          defaults to 0.5
        :param include_bond_order: If bond orders should be included as edge
                                   attributes, defaults to False
        :return: Returns MolGraph
        """
        # scipy is not a dependency, this import is only done if needed
        # so that molgraph can be used without scipy
        # and to suppress the warning
        
        import scipy

        if not len(atom_types) == np.shape(sp_arr)[0] == np.shape(sp_arr)[1]:
            raise ValueError(
                "atom_types and matrix have to have the same length"
            )
        new_mol_graph = cls()
        with warnings.catch_warnings():
            warnings.simplefilter(
                action="ignore", category=scipy.sparse.SparseEfficiencyWarning
            )

            # catching SparseEfficiencyWarning:
            # Comparing a sparse matrix with a scalar greater than zero
            # using < is inefficient, try using >= instead.
            gt_thresh = (threshold > sp_arr) * sp_arr

        if include_bond_order is True:
            new_graph = new_mol_graph._graph.from_numpy_array(
                gt_thresh, edge_attr="bond_order"
            )
        else:
            new_graph = new_mol_graph._graph.from_numpy_array(gt_thresh)

        new_graph.set_node_attributes(
            {i: PERIODIC_TABLE[atom_type] for i, atom_type in enumerate(atom_types)},
            name="atom_type",
        )

        new_mol_graph._graph = new_graph

        return new_mol_graph

    @classmethod
    def from_geometry_and_bond_order_matrix(
        cls: type[Self],
        geo: Geometry,
        matrix: np.ndarray,
        threshold: float = 0.5,
        include_bond_order: bool = False,
    ) -> Self:
        """
        Creates a graph of a molecule from a Geometry and a bond order matrix.

        :param geo: Geometry
        :param matrix: Bond order matrix
        :param threshold: Threshold for bonds to be included as edges,
            defaults to 0.5
        :param include_bond_order: If bond orders should be included as edge
            attributes, defaults to False
        :return: Graph of Molecule
        """
        new_mol_graph = cls.from_atom_types_and_bond_order_matrix(
            geo.atom_types,
            matrix,
            threshold=threshold,
            include_bond_order=include_bond_order,
        )
        return new_mol_graph

    @classmethod
    def from_geometry(
        cls: type[Self],
        geo: Geometry,
        switching_function: Callable = BondsFromDistance(),
    ) -> Self:
        """
        Creates a graph of a molecule from a Geometry and a switching Function.
        Uses the Default switching function if none are given.

        :param geo: Geometry
        :param switching_function: Function to determine if two atoms are
            connected
        :return: graph of molecule
        """

        connectivity_matrix = switching_function.array(
            geo.coords, geo.atom_types
        )
        return cls.from_geometry_and_bond_order_matrix(
            geo,
            connectivity_matrix,
        )

    def __eq__(self, other: object) -> bool:
        if isinstance(other, MolGraph):
            return self.is_isomorphic(other)
        return NotImplemented

    def __hash__(self) -> int:
        return self.color_refine_hash()

    def color_refine_hash(self) -> int:
        """TODO"""
        color_dict = color_refine_mg(self, )
        return hash(tuple(sorted(Counter(color_dict.values()).items())))

    def bfs_layers(self, sources: Iterable[AtomId] | AtomId
                   ) -> Iterator[list[AtomId]]:
        """
        Generates layers of the graph starting from the source atoms.
        Each layer contains all nodes that are at the same distance from the
        sources.
        The first layer contains the sources.
        
        :param sources: Sources to start from
        """
        if sources in self.atoms:
            sources = [sources]

        current_layer = list(sources)
        visited = set(sources)

        if any(source not in self.atoms for source in current_layer):
            raise ValueError("Source atom not in graph")

        # this is basically BFS, except that the current layer only stores
        # the nodes at same distance from sources at each iteration
        while current_layer:
            yield current_layer
            next_layer = []
            for node in current_layer:
                for child in self.bonded_to(node):
                    if child not in visited:
                        visited.add(child)
                        next_layer.append(child)
            current_layer = next_layer

    def get_subgraph_isomorphic_mappings(
        self, other: Self
    ) -> Iterator[dict[int, int]]:
        """Subgraph isomorphic mappings from "other" onto "self".

        Generates all node-iduced subgraph isomorphic mappings.
        All atoms of "other" have to be present in "self".
        The bonds of "other" have to be the subset of the bonds of "self"
        relating to the nodes of "other".

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        """
        return vf2pp_all_isomorphisms(
            self,
            other,
            color_refine=False,
            stereo=False,
            stereo_change=False,
            subgraph=True,
        )

    def get_isomorphic_mappings(self, other: Self) -> Iterator[dict[int, int]]:
        """Isomorphic mappings between "self" and "other".

        Generates all isomorphic mappings between "other" and "self".
        All atoms and bonds have to be present in both graphs.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """

        return vf2pp_all_isomorphisms(
            self,
            other,
            color_refine=False, # TODO: implement color refinement
            stereo=False,
            stereo_change=False,
            subgraph=False,
        )

    def is_isomorphic(self, other: Self) -> bool:
        """
        Checks if the graph is isomorphic to another graph.

        :param other: other graph
        :return: True if isomorphic
        """
        return any(self.get_isomorphic_mappings(other))


class BondChange(Enum):
    FORMED = 1
    FLEETING = 0
    BROKEN = -1

    def __repr__(self):
        return self.name


class CondensedReactionGraph(MolGraph):
    """
    Graph representing a reaction. Atoms are nodes and (potentially changing)
    bonds are edges. Every node has to have an attribute "atom_type" of type
    Element. Edges can have an attribute "reaction" of type BondChange.
    This is used to represent the change in connectivity during the reaction.

    Two graphs are equal, iff. they are isomporhic and of the same type.
    """
    __slots__ = tuple()
    _atom_attrs: dict[AtomId, dict[str, Any]]
    _neighbors: dict[AtomId, set[AtomId]]
    _bond_attrs: dict[Bond, dict[str, Any]]

    def add_bond(self, atom1: int, atom2: int, **attr):
        """
        Adds a bond between atom1 and atom2.

        :param atom1: id of atom1
        :param atom2:   id of atom2
        """
        if "reaction" in attr and not isinstance(
            attr.get("reaction"), BondChange
        ):
            raise TypeError("reaction bond has to have reaction attribute")
        super().add_bond(atom1, atom2, **attr)

    def set_bond_attribute(
        self, atom1: int, atom2: int, attr: str, value: Any
    ):
        """
        sets the Attribute of the bond between Atom1 and Atom2.

        :param atom1: Atom1
        :param atom2: Atom2
        :param attr: Attribute
        :param value: Value
        """
        if attr == "reaction" and not isinstance(value, BondChange):
            raise ValueError("reaction bond has to have reaction attribute")
        super().set_bond_attribute(atom1, atom2, attr, value)

    def add_formed_bond(self, atom1: int, atom2: int, **attr):
        """
        Adds a bond between atom1 and atom2 with reaction attribute
        set to FORMED.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        if atom1 in self._atom_attrs and atom2 in self._atom_attrs:
            self.add_bond(atom1, atom2, reaction=BondChange.FORMED, **attr)
        else:
            raise ValueError("Atoms have to be in the graph")

    def add_broken_bond(self, atom1: int, atom2: int, **attr):
        """
        Adds a bond between atom1 and atom2 with reaction attribute
        set to BROKEN.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        if atom1 in self._atom_attrs and atom2 in self._atom_attrs:
            self.add_bond(atom1, atom2, reaction=BondChange.BROKEN, **attr)
        else:
            raise ValueError("Atoms have to be in the graph")
            

    def get_formed_bonds(self) -> set[Bond[int, int]]:
        """
        Returns all bonds that are formed during the reaction

        :return: formed bonds
        """
        return {
            bond
            for bond in self.bonds
            if self.get_bond_attribute(*bond, "reaction") == BondChange.FORMED
        }

    def get_broken_bonds(self) -> set[Bond[int, int]]:
        """
        Returns all bonds that are broken during the reaction

        :return: broken bonds
        """
        return {
            bond
            for bond in self.bonds
            if self.get_bond_attribute(*bond, "reaction") == BondChange.BROKEN
        }

    def active_atoms(self, additional_layer: int = 0) -> set[int]:
        """
        Atoms involved in the reaction with additional layers of atoms
        in the neighborhood.

        :param additional_layer: Number of additional layers of atoms to
                                 include, defaults to 0
        :return: Atoms involved in the reaction
        """
        active_atoms: set[int] = set()
        for bond in self.get_formed_bonds() | self.get_broken_bonds():
            active_atoms.update(bond)
        for _ in range(additional_layer):
            for atom in active_atoms.copy():
                active_atoms.update(self._neighbors[atom])
        return active_atoms

    def connectivity_matrix(
        self,
    ) -> np.ndarray:
        """
        Returns a connectivity matrix of the graph. Order is the same
        as in self.atoms
        Formed bonds and broken bonds are represented as 0.5.

        :return: Connectivity matrix
        """

        matrix = np.array(super().connectivity_matrix(), dtype=float)
        atoms = tuple(self.atoms)
        for bond in self.get_formed_bonds():
            a1, a2 = bond
            index_atom1 = atoms.index(a1)
            index_atom2 = atoms.index(a2)

            matrix[index_atom1][index_atom2] = 0.5
            matrix[index_atom2][index_atom1] = 0.5

        for bond in self.get_broken_bonds():
            a1, a2 = bond
            index_atom1 = atoms.index(a1)
            index_atom2 = atoms.index(a2)

            matrix[index_atom1][index_atom2] = 0.5
            matrix[index_atom2][index_atom1] = 0.5
        return matrix

    def _to_rdmol(
        self, generate_bond_orders=False, charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        mol, idx_map_num_dict = _mol_graph_to_rdmol(graph=self, generate_bond_orders=False, charge=0)
        _set_crg_bond_orders(graph=self, mol=mol, idx_map_num_dict=idx_map_num_dict)
        return mol, idx_map_num_dict

    def to_rdmol(self) -> Chem.rdchem.RWMol:
        raise NotImplementedError(
            "Rdkit is not able to represent "
            "reactions as condensed reaction graphs."
        )

    def reactant(self, keep_attributes=True) -> MolGraph:
        """Reactant of the reaction

        Creates the reactant of the reaction.
        Formed bonds are not present in the reactant.

        :param keep_attributes: attributes on atoms and bonds to be kept,
                                defaults to True
        :return: Reactant of the reaction
        """
        product = MolGraph()
        for atom in self.atoms:
            if keep_attributes is True:
                attrs = self._atom_attrs[atom]
            else:
                attrs = {
                    "atom_type": self._atom_attrs[atom]["atom_type"]
                }
            product.add_atom(atom, **attrs)
        for bond in self.bonds:
            bond_reaction = self._bond_attrs[bond].get("reaction", None)
            if (
                bond_reaction is None or bond_reaction == BondChange.BROKEN
            ):
                if keep_attributes is True:
                    attrs = self._bond_attrs[bond].copy()
                    attrs.pop("reaction", None)
                else:
                    attrs = {}
                product.add_bond(*bond, **attrs)
        return product

    def product(self, keep_attributes=True) -> MolGraph:
        """Product of the reaction

        Creates the product of the reaction.
        Broken bonds are not present in the product.

        :param keep_attributes: attributes on atoms and bonds to be kept,
                                defaults to True
        :return: Product of the reaction
        """
        product = MolGraph()
        for atom in self.atoms:
            if keep_attributes is True:
                attrs = self._atom_attrs[atom]
            else:
                attrs = {
                    "atom_type": self._atom_attrs[atom]["atom_type"]
                }
            product.add_atom(atom, **attrs)
        for bond in self.bonds:
            bond_reaction = self._bond_attrs[bond].get("reaction", None)
            if (
                bond_reaction is None or bond_reaction == BondChange.FORMED
            ):
                if keep_attributes is True:
                    attrs = self._bond_attrs[bond].copy()
                    attrs.pop("reaction", None)
                else:
                    attrs = {}
                product.add_bond(*bond, **attrs)
        return product

    def reverse_reaction(self) -> Self:
        """Creates the reaction in the opposite direction.

        Broken bonds are turned into formed bonds and the other way around.

        :return: Reversed reaction
        """
        rev_reac = self.copy()
        for bond in self.bonds:
            bond_reaction = self._bond_attrs[bond].get("reaction", None)
            if bond_reaction == BondChange.FORMED:
                rev_reac.add_broken_bond(*bond)
            elif bond_reaction == BondChange.BROKEN:
                rev_reac.add_formed_bond(*bond)

        return rev_reac

    @classmethod
    def from_reactant_and_product_graph(
        cls: type[Self], reactant_graph: MolGraph, product_graph: MolGraph
    ) -> Self:
        """Creates a CondensedReactionGraph from reactant and product MolGraphs

        CondensedReactionGraph  is constructed from bond changes from reactant
        to the product. The atoms order and atom types of the reactant and
        product have to be the same.

        :param reactant_graph: reactant of the reaction
        :param product_graph: product of the reaction
        :return: CondensedReactionGraph
        """

        if set(reactant_graph.atoms) != set(product_graph.atoms):
            raise ValueError("reactant and product have different atoms1")
        for atom in reactant_graph.atoms:
            if reactant_graph.get_atom_attribute(
                atom, "atom_type"
            ) != product_graph.get_atom_attribute(atom, "atom_type"):
                raise ValueError(
                    "reactant and product have different atom types"
                )

        crg = cls()

        atoms = {*reactant_graph.atoms, *product_graph.atoms}
        for atom in atoms:
            crg.add_atom(
                atom,
                atom_type=reactant_graph.get_atom_attribute(atom, "atom_type"),
            )

        bonds = {
            *[tuple(sorted(bond)) for bond in reactant_graph.bonds],
            *[tuple(sorted(bond)) for bond in product_graph.bonds],
        }

        for bond in bonds:
            if reactant_graph.has_bond(*bond) and product_graph.has_bond(
                *bond
            ):
                crg.add_bond(*bond)
            elif reactant_graph.has_bond(*bond):
                crg.add_bond(*bond, reaction=BondChange.BROKEN)
            elif product_graph.has_bond(*bond):
                crg.add_bond(*bond, reaction=BondChange.FORMED)
        return crg

    @classmethod
    def from_reactant_and_product_geometry(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        switching_function: Callable = BondsFromDistance(),
    ) -> Self:
        """Creates a CondensedReactionGraph from reactant
        and product Geometries.


        CondensedReactionGraph  is constructed from bond changes from reactant
        to the product. The atoms order and atom types of the reactant and
        product have to be the same. The switching function is used to
        determine the connectivity of the atoms.

        :param reactant_geo: geometry of the reactant
        :param product_geo: geometry of the product
        :param switching_function: function to define the connectivity
                                   from geometry,
                                   defaults to StepSwitchingFunction()
        :return: CondensedReactionGraph
        """

        reactant = MolGraph.from_geometry(reactant_geo, switching_function)
        product = MolGraph.from_geometry(product_geo, switching_function)
        return cls.from_reactant_and_product_graph(reactant, product)

    def get_isomorphic_mappings(self, other: Self) -> Iterator[dict[int, int]]:
        """Isomorphic mappings between "self" and "other".

        Generates all isomorphic mappings between "other" and "self".
        All atoms and bonds have to be present in both graphs.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """

        return vf2pp_all_isomorphisms(
            self,
            other,
            color_refine=False, # TODO: implement color refinement
            stereo=False,
            stereo_change=False,
            subgraph=False,
            # labels=["reaction"],
        )

    def apply_reaction(
        self,
        reactant: MolGraph,
        mapping: Mapping[AtomId, AtomId],
    ) -> Self:
        """
        Applies a reaction to the graph and returns the resulting graph.
        The reactants of the CRG have to be a subgraph of the reactant.
        Mappings from crg atoms to the reactant can be provided.
        Atom numberig from the reactant is kept in the resulting graph.

        :param reaction: Reaction to apply
        :param mapping: Mappings from reactant atoms to crg atoms
        :return: Resulting graph
        """
        crg = self.__class__(reactant)

        for a1, a2 in self.get_formed_bonds():
            crg.add_formed_bond(mapping[a1], mapping[a2])
        for a1, a2 in self.get_broken_bonds():
            crg.remove_bond(mapping[a1], mapping[a2])
            crg.add_broken_bond(mapping[a1], mapping[a2])

        return crg


class StereoMolGraph(MolGraph):
    """
    :class:`MolGraph` with the ability to store stereochemistry information
    for atoms and bonds.

    Two graphs compare equal, if they are isomorphic and have the same
    stereochemistry.
    """
    __slots__ = ("_atom_stereo", "_bond_stereo")
    _atom_stereo: dict[int, AtomStereo]
    _bond_stereo: dict[Bond, PlanarBond]

    def __init__(self, mol_graph: Optional[MolGraph] = None):
        super().__init__(mol_graph)
        if mol_graph and isinstance(mol_graph, StereoMolGraph):
            self._atom_stereo = deepcopy(mol_graph._atom_stereo)
            self._bond_stereo = deepcopy(mol_graph._bond_stereo)
        else:
            self._atom_stereo = {}
            self._bond_stereo = {}

    @property
    def stereo(self) -> Mapping[AtomId | Bond, AtomStereo | BondStereo]:
        return MappingProxyType(self._atom_stereo | self._bond_stereo)

    @property
    def atom_stereo(self) -> Mapping[AtomId, AtomStereo]:
        return MappingProxyType(self._atom_stereo)

    @property
    def bond_stereo(self) -> Mapping[Bond, BondStereo]:
        return MappingProxyType(self._bond_stereo)

    def get_atom_stereo(
        self, atom: AtomId
    ) -> Optional[AtomStereo]:
        """Returns the stereo information of the atom if it exists else None.
        Raises a ValueError if the atom is not in the graph.

        :param atom: atom
        :param default: Default value if no stereo information is found,
                        defaults to None
        :return: Stereo information of atom
        """
        if atom in self._atom_attrs:
            if s := self._atom_stereo.get(atom, None):
                return s
            else:
                return None
                #return NoStereo(atoms=(atom, *list(self.bonded_to(atom))))
        else:
            raise ValueError(f"Atom {atom} is not in the graph")

    def set_atom_stereo(self, atom: AtomId, atom_stereo: AtomStereo):
        """Adds stereo information to the graph

        :param atom: Atoms to be used for chiral information
        :param stereo: Chiral information
        """
        if atom in self._atom_attrs:
            assert atom in atom_stereo.atoms
            self._atom_stereo[atom] = atom_stereo
        else:
            raise ValueError(f"Atom {atom} is not in the graph")

    def delete_atom_stereo(self, atom: AtomId):
        """Deletes stereo information from the graph

        :param atom: Atom to be used for stereo information
        """
        del self._atom_stereo[atom]

    def get_bond_stereo(
        self, bond: Iterable[int]
    ) -> Optional[BondStereo]:
        """Gets the stereo information of the bond or None
        if it does not exist.
        Raises a ValueError if the bond is not in the graph.

        :param bond: Bond
        :return: stereo information of bond
        """
        bond = Bond(bond)
        bond_stereo = self._bond_stereo.get(Bond(bond), None)
        if bond_stereo:
            return bond_stereo
        elif bond in self._bond_attrs:
            return None
        else:
            raise ValueError(f"Bond {bond} is not in the graph")

    def set_bond_stereo(
        self, bond: Iterable[int], bond_stereo: BondStereo
    ):
        """Stets the stereo information of the bond

        :param bond: Bond
        :param bond_stereo: Stereo information of the bond
        """
        bond = Bond(bond)
        if bond in self._bond_attrs:
            self._bond_stereo[Bond(bond)] = bond_stereo
        else:
            raise ValueError(f"Bond {bond} is not in the graph")

    def delete_bond_stereo(self, bond: Iterable[int]):
        """Deletes the stereo information of the bond

        :param bond: Bond
        """
        del self._bond_stereo[Bond(bond)]

    def delete_stereo(self, atom_or_bond: AtomId | Iterable[AtomId]):
        """Deletes the stereo information of the atom or bond

        :param atom_or_bond: Atom or Bond
        """
        if isinstance(atom_or_bond, int):
            self.delete_atom_stereo(atom_or_bond)
        elif isinstance(atom_or_bond, Iterable):
            self.delete_bond_stereo(atom_or_bond)
        else:
            raise TypeError("atom_or_bond and stereo have the wrong type")

    def remove_atom(self, atom: int):
        """Removes an atom from the graph and deletes all chiral information
        associated with it

        :param atom: Atom
        """
        for a, atom_stereo in self._atom_stereo.copy().items():
            if atom in atom_stereo.atoms:
                self.delete_atom_stereo(a)

        for bond, bond_stereo in self._bond_stereo.copy().items():
            if atom in bond_stereo.atoms:
                self.delete_bond_stereo(bond)
        super().remove_atom(atom)

    def copy(self) -> Self:
        """
        :return: returns a copy of self
        """
        new_graph = super().copy()
        new_graph._atom_stereo = deepcopy(self._atom_stereo)
        new_graph._bond_stereo = deepcopy(self._bond_stereo)
        return new_graph

    def relabel_atoms(
        self, mapping: dict[int, int], copy: bool = True
    ) -> Self:
        """
        Relabels the atoms of the graph and the chiral information accordingly

        :param mapping: Mapping of old atom ids to new atom ids
        :param copy: If the graph should be copied before relabeling,
                     defaults to True
        :return: Returns the relabeled graph
        """
        new_atom_stereo_dict = self._atom_stereo.__class__()
        new_bond_stereo_dict = self._bond_stereo.__class__()

        for central_atom, stereo in self._atom_stereo.items():
            new_central_atom = mapping.get(central_atom, central_atom)
            new_atom_stereo_atoms = tuple(
                mapping.get(atom, atom) for atom in stereo.atoms
            )
            new_atom_stereo = stereo.__class__(
                new_atom_stereo_atoms, stereo.parity
            )
            new_atom_stereo_dict[new_central_atom] = new_atom_stereo

        for bond, bond_stereo in self._bond_stereo.items():
            new_bond = tuple(mapping.get(atom, atom) for atom in bond)
            new_bond_stereo_atoms = tuple(
                mapping.get(atom, atom) for atom in bond_stereo.atoms
            )
            new_bond_stereo = bond_stereo.__class__(
                new_bond_stereo_atoms, bond_stereo.parity
            )
            new_bond_stereo_dict[frozenset(new_bond)] = new_bond_stereo

        if copy is True:
            graph = super().relabel_atoms(mapping, copy=True)
            graph._atom_stereo = new_atom_stereo_dict
            graph._bond_stereo = new_bond_stereo_dict
            return graph

        elif copy is False:
            super().relabel_atoms(mapping, copy=False)
            self._atom_stereo = new_atom_stereo_dict
            self._bond_stereo = new_bond_stereo_dict
            return self

    def subgraph(self, atoms: Iterable[int]) -> Self:
        """Returns a subgraph of the graph with the given atoms and the chiral
        information accordingly

        :param atoms: Atoms to be used for the subgraph
        :return: Subgraph
        """
        new_graph = super().subgraph(atoms)

        for central_atom, atoms_atom_stereo in self._atom_stereo.items():
            atoms_set = set((*atoms_atom_stereo.atoms, central_atom))
            if all(atom in atoms for atom in atoms_set):
                new_graph.set_atom_stereo(central_atom, atoms_atom_stereo)

        for bond, bond_stereo in self._bond_stereo.items():
            if all(atom in atoms for atom in bond_stereo.atoms):
                new_graph.set_bond_stereo(bond, bond_stereo)
        return new_graph

    def enantiomer(self) -> Self:
        """
        Creates the enantiomer of the StereoMolGraph by inversion of all atom
        stereocenters. The result can be identical to the molecule itself if
        no enantiomer exists.

        :return: Enantiomer
        """
        enantiomer = self.copy()
        for atom in self.atoms:
            if stereo := self.get_atom_stereo(atom):
                enantiomer.set_atom_stereo(atom, stereo.invert())
        return enantiomer

    def _to_rdmol(
        self, generate_bond_orders=False, charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        """
        Creates a RDKit mol object using the connectivity of the mol graph.
        Stereochemistry is added to the mol object.

        :return: RDKit molecule
        """
        return _stereo_mol_graph_to_rdmol(self,
                                          generate_bond_orders=generate_bond_orders,
                                          charge=charge)
        

    @classmethod
    def from_rdmol(cls, rdmol, use_atom_map_number=False) -> Self:
        """
        Creates a StereoMolGraph from an RDKit Mol object.
        All hydrogens have to be explicit.
        Stereo information is conserved for tetrahedral atoms and
        double bonds.

        :param rdmol: RDKit Mol object
        :param use_atom_map_number: If the atom map number should be used
                                    instead of the atom index, Default: False
        :return: StereoMolGraph
        """
        return stereo_mol_graph_from_rdmol(cls, rdmol, use_atom_map_number=use_atom_map_number)

    def _set_atom_stereo_from_geometry(self, geo: Geometry):
        for atom in range(geo.n_atoms):
            first_neighbors = self.bonded_to(atom)

            # extends the first layer of neighbors to the second layer
            # (if planar)
            # needed to find double bonds
            if len(first_neighbors) < 3:
                pass
            elif len(first_neighbors) == 3 and are_planar(
                *(
                    geo.coords[i]
                    for i in (
                        list(first_neighbors)
                        + [
                            atom,
                        ]
                    )
                )
            ):
                for neighbor in first_neighbors:
                    next_layer = set(self.bonded_to(neighbor))

                    for i in first_neighbors:
                        next_layer.add(i)
                    next_layer -= set((neighbor, atom))

                    if len(next_layer) != 4:
                        continue

                    elif are_planar(*(geo.coords[i] for i in next_layer)):
                        bonded_to_atom = (
                            outer_atom
                            for outer_atom in next_layer
                            if self.has_bond(atom, outer_atom)
                        )
                        bonded_to_neighbor = (
                            outer_atom
                            for outer_atom in next_layer
                            if self.has_bond(neighbor, outer_atom)
                        )
                        atom_ids = (
                            *bonded_to_atom,
                            atom,
                            neighbor,
                            *bonded_to_neighbor,
                        )
                        double_bond = PlanarBond.from_coords(
                            atom_ids, *tuple(geo.coords[i] for i in atom_ids)
                        )
                        self.set_bond_stereo((atom, neighbor), double_bond)

            elif len(first_neighbors) == 3:
                pass
            else:
                first_neighbors_coords = tuple(
                    geo.coords[i] for i in first_neighbors
                )
                if are_planar(*first_neighbors_coords):
                    pass

                elif len(first_neighbors) == 4:
                    atoms_atom_stereo = Tetrahedral.from_coords(
                        (atom, *first_neighbors), None, *first_neighbors_coords
                    )
                    self.set_atom_stereo(atom, atoms_atom_stereo)

                elif len(first_neighbors) == 5:
                    atoms_atom_stereo = TrigonalBipyramidal.from_coords(
                        (atom, *first_neighbors), None, *first_neighbors_coords
                    )
                    self.set_atom_stereo(atom, atoms_atom_stereo)

    @classmethod
    def from_composed_molgraphs(cls, mol_graphs: Iterable[Self]) -> Self:
        """Creates a MolGraph object from a list of MolGraph objects.
        
        Duplicate nodes or edges are overwritten, such that the resulting
        graph only contains one node or edge with that name. Duplicate
        attributes of duplicate nodes, edges and the stereochemistry are also
        overwritten in order of iteration.

        :param mol_graphs: list of MolGraph objects
        :return: Returns MolGraph
        """

        graph = cls(super().from_composed_molgraphs(mol_graphs))
        for mol_graph in mol_graphs:
            graph._atom_stereo.update(mol_graph._atom_stereo)
            graph._bond_stereo.update(mol_graph._bond_stereo)
        return graph

    @classmethod
    def from_geometry_and_bond_order_matrix(
        cls: type[Self],
        geo: Geometry,
        matrix: np.ndarray,
        threshold: float = 0.5,
        include_bond_order: bool = False,
    ) -> Self:
        """
        Creates a CiralMolGraph object from a Geometry and a bond order matrix

        :param geo: Geometry
        :param matrix: Bond order matrix
        :param threshold: Threshold for bonds to be included as edges,
                          defaults to 0.5
        :param include_bond_order: If bond orders should be included as edge
                                    attributes, defaults to False
        :return: Returns MolGraph
        """
        mol_graph = super().from_geometry_and_bond_order_matrix(
            geo,
            matrix=matrix,
            threshold=threshold,
            include_bond_order=include_bond_order,
        )
        graph = cls(mol_graph)
        graph._set_atom_stereo_from_geometry(geo)
        return graph

    def get_isomorphic_mappings(
        self, other: Self, stereo=True
    ) -> Iterator[dict[int, int]]:
        """Isomorphic mappings between "self" and "other".

        Generates all isomorphic mappings between "other" and "self".
        All atoms and bonds have to be present in both graphs.
        The Stereochemistry is preserved in the mappings.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """

        return vf2pp_all_isomorphisms(
            self,
            other,
            color_refine=False, #TODO: implement color refinement
            stereo=stereo,
            stereo_change=False,
            subgraph=False,
        )
        

    def get_subgraph_isomorphic_mappings(
        self, other: Self, stereo: bool = True
    ) -> Iterator[dict[int, int]]:
        """Subgraph isomorphic mappings from "self" onto "other".
        Other can be of equal size or larger than "self".
        Generates all node-iduced subgraph isomorphic mappings.
        All atoms of "self" have to be present in "other".
        The bonds of "self" have to be the subset of the bonds of "other"
        relating to the nodes of "self".
        The Stereochemistry is preserved in the mappings.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """
        return vf2pp_all_isomorphisms(
            self,
            other,
            color_refine=False,
            stereo=stereo,
            stereo_change=False,
            subgraph=True,
        )



    def is_stereo_valid(self) -> bool:
        """
        Checks if the bonds required to have the defined stereochemistry
        are present in the graph.

        :return: True if the stereochemistry is valid
        """
        for atom, stereo in self._atom_stereo.items():
            for neighbor in stereo.atoms[1:]:
                if not self.has_bond(atom, neighbor):
                    return False
        for bond, stereo in self._bond_stereo.items():
            if not self.has_bond(*bond):
                return False
            if {stereo.atoms[2], stereo.atoms[3]} != set(bond):
                return False
            if not self.has_bond(stereo.atoms[0], stereo.atoms[2]):
                return False
            if not self.has_bond(stereo.atoms[1], stereo.atoms[2]):
                return False
            if not self.has_bond(stereo.atoms[4], stereo.atoms[3]):
                return False
            if not self.has_bond(stereo.atoms[5], stereo.atoms[3]):
                return False
        return True


class StereoChange(Enum):
    BROKEN = "broken"
    FLEETING = "fleeting"
    FORMED = "formed"

    def __repr__(self):
        return self.name


class StereoChangeDict(dict[StereoChange, Stereo]):
    def __missing__(self, key: StereoChange):
        if key in StereoChange:
            return None
        else:
            raise KeyError(f"{key} not in {self.__class__.__name__}")
    

class StereoCondensedReactionGraph(StereoMolGraph, CondensedReactionGraph):
    """
    :class:`CondenedReactionGraph` with the ability to store stereochemistry
    information for atoms and (potentially changing) bonds.
    """

    __slots__ = ("_atom_stereo_change", "_bond_stereo_change")
    _atom_stereo_change: Mapping[AtomId, Mapping[StereoChange, AtomStereo]]
    _bond_stereo_change: Mapping[Bond, Mapping[StereoChange, BondStereo]]

    def __init__(self, mol_graph: Optional[MolGraph] = None):
        super().__init__(mol_graph)
        if mol_graph and isinstance(mol_graph, StereoCondensedReactionGraph):
            self._atom_stereo_change = mol_graph._atom_stereo_change.copy()
            self._bond_stereo_change = mol_graph._bond_stereo_change.copy()
        else:
            self._atom_stereo_change = defaultdict(StereoChangeDict)
            self._bond_stereo_change = defaultdict(StereoChangeDict)

    @property
    def atom_stereo_changes(self) -> Mapping[AtomId, StereoChangeDict]:
        return MappingProxyType(self._atom_stereo_change)

    @property
    def bond_stereo_changes(self) -> Mapping[Bond, StereoChangeDict]:
        return MappingProxyType(self._bond_stereo_change)

    @property
    def stereo_changes(self) -> Mapping[AtomId | Bond, StereoChangeDict]:
        return MappingProxyType(self._atom_stereo_change
                                | self._bond_stereo_change)

    def get_atom_stereo_change(
        self, atom: int
    ) -> Mapping[StereoChange, AtomStereo]:
        
        if atom in self._atom_attrs:
            if atom in self._atom_stereo_change:
                return MappingProxyType(self._atom_stereo_change[atom])
            else:
                return None
        else:
            raise ValueError(f"Atom {atom} not in graph")

    def get_bond_stereo_change(
        self, bond: Iterable[int]
        ) -> Mapping[StereoChange, BondStereo]:
        
        bond = Bond(bond)
        if bond in self._bond_attrs:
            if bond in self._bond_stereo_change:
                return MappingProxyType(self._bond_stereo_change[bond])
            else:
                return None
        else:
            raise ValueError(f"Bond {bond} not in graph")
        
    def set_atom_stereo_change(
        self,
        atom: AtomId,
        *,
        broken: Optional[AtomStereo] = None,
        fleeting: Optional[AtomStereo] = None,
        formed: Optional[AtomStereo] = None,
    ):
        if atom not in self._atom_attrs:
            raise ValueError(f"Atom {atom} not in graph")
        for stereo_change, atom_stereo in {
            StereoChange.BROKEN: broken,
            StereoChange.FLEETING: fleeting,
            StereoChange.FORMED: formed,
        }.items():
            if atom_stereo:
                self._atom_stereo_change[atom][stereo_change] = atom_stereo

    def set_bond_stereo_change(
        self,
        bond: Iterable[AtomId],
        broken: Optional[BondStereo] = None,
        fleeting: Optional[BondStereo] = None,
        formed: Optional[BondStereo] = None,

    ):
        bond = Bond(bond)
        if bond not in self._bond_attrs:
            raise ValueError(f"Bond {bond} not in graph")
        bond = Bond(bond)
        for stereo_change, bond_stereo in {
            StereoChange.BROKEN: broken,
            StereoChange.FORMED: formed,
            StereoChange.FLEETING: fleeting
        }.items():
            if bond_stereo:
                self._bond_stereo_change[bond][stereo_change] = bond_stereo

    def delete_atom_stereo_change(
        self, atom: AtomId, stereo_change: Optional[StereoChange] = None
    ):
        if stereo_change is None:
            del self._atom_stereo_change[atom]
        else:
            del self._atom_stereo_change[atom][stereo_change]

    def delete_bond_stereo_change(
        self, bond: Iterable[AtomId],
        stereo_change: Optional[StereoChange] = None
    ):
        bond = Bond(bond)
        if stereo_change is None:
            del self._bond_stereo_change[bond]
        else:
            del self._bond_stereo_change[bond][stereo_change]

    def active_atoms(self, additional_layer: int = 0) -> set[AtomId]:
        """
        Atoms involved in the reaction with additional layers of atoms
        in the neighborhood.

        :param additional_layer: Number of additional layers of atoms to
                                 include, defaults to 0
        :return: Atoms involved in the reaction
        """
        active_atoms: set[int] = set()
        for bond in self.get_formed_bonds() | self.get_broken_bonds():
            active_atoms.update(bond)
        for atom_or_bond, stereo_change in self.stereo_changes.items():
            for change, stereo in stereo_change.items():
                for atom in stereo.atoms:
                    active_atoms.update(atom)
        for _ in range(additional_layer):
            for atom in active_atoms.copy():
                active_atoms.update(self.bonded_to(atom))
        return active_atoms

    def copy(self) -> Self:
        """
        :return: returns a copy of self
        """
        new_graph = super().copy()
        new_graph._atom_stereo_change = deepcopy(self._atom_stereo_change)
        new_graph._bond_stereo_change = deepcopy(self._bond_stereo_change)
        return new_graph

    def relabel_atoms(
        self, mapping: dict[AtomId, AtomId], copy: bool = True
    ) -> Self:
        """
        Relabels the atoms of the graph and the chiral information accordingly

        :param mapping: Mapping of old atom ids to new atom ids
        :param copy: If the graph should be copied before relabeling,
                     defaults to True
        :return: Returns the relabeled graph or None if copy is False
        """
        relabeled_scrg = self.__class__(
            super().relabel_atoms(mapping, copy=copy)
        )

        atom_stereo_change = defaultdict(StereoChangeDict)
        
        for atom, stereo_change_dict in self._atom_stereo_change.items():
            for stereo_change, atom_stereo in stereo_change_dict.items():
                new_stereo = atom_stereo.__class__(
                    atoms = tuple(mapping.get(atom, atom)
                                  for atom in atom_stereo.atoms),
                    parity = atom_stereo.parity,
                )
                atom_stereo_change[mapping[atom]][stereo_change] = new_stereo

        bond_stereo_change = defaultdict(StereoChangeDict)
        
        for bond, stereo_change_dict in self._bond_stereo_change.items():
            for stereo_change, bond_stereo in stereo_change_dict.items():
                new_bond = Bond((mapping[bond[0]], mapping[bond[1]]))
                new_stereo = bond_stereo.__class__(
                    atoms = tuple(mapping.get(atom, atom)
                                  for atom in bond_stereo.atoms),
                    parity = bond_stereo.parity,
                )
            bond_stereo_change[new_bond][stereo_change] = new_stereo
            
        if copy is True:
            relabeled_scrg._atom_stereo_change = atom_stereo_change
            relabeled_scrg._bond_stereo_change = bond_stereo_change
        else:
            self._atom_stereo_change = atom_stereo_change
            self._bond_stereo_change = bond_stereo_change

        return relabeled_scrg

    def reactant(self, keep_attributes=True) -> StereoMolGraph:
        """
        Returns the reactant of the reaction

        :param keep_attributes: If attributes should be kept , defaults to True
        :return: reactant
        """

        reactant = StereoMolGraph(
            super().reactant(keep_attributes=keep_attributes)
        )
        reactant._atom_stereo = deepcopy(self._atom_stereo)
        reactant._bond_stereo = deepcopy(self._bond_stereo)

        for atom, change_dict in self._atom_stereo_change.items():
            if stereo := change_dict[StereoChange.BROKEN]:
                reactant._atom_stereo[atom] = stereo
                

        for bond, change_dict in self._bond_stereo_change.items():
            if stereo := change_dict[StereoChange.BROKEN]:
                reactant._bond_stereo[bond] = stereo
                #reactant.set_bond_stereo(bond, stereo)

        return reactant

    def product(self, keep_attributes=True) -> StereoMolGraph:
        """
        Returns the product of the reaction

        :param keep_attributes: If attributes should be kept, defaults to True
        :return: product
        """
        product = StereoMolGraph(
            super().product(keep_attributes=keep_attributes)
        )
        product._atom_stereo = deepcopy(self._atom_stereo)
        product._bond_stereo = deepcopy(self._bond_stereo)

        for atom, change_dict in self._atom_stereo_change.items():
            if stereo := change_dict[StereoChange.FORMED]:
                product.set_atom_stereo(atom, stereo)

        for bond, change_dict in self._bond_stereo_change.items():
            if stereo := change_dict[StereoChange.FORMED]:
                product.set_bond_stereo(bond, stereo)

        return product

    def reverse_reaction(self) -> Self:
        """Creates the reaction in the opposite direction.

        Broken bonds and stereochemistry changes are turned into formed
        and the other way around.

        :return: Reversed reaction
        """
        rev_reac = super().reverse_reaction()
        for atom, change_dict in rev_reac._atom_stereo_change.items():
            new_change_dict = {
                stereo_change.value: atom_stereo
                for stereo_change, atom_stereo in change_dict.items()
            }
            if formed_stereo := change_dict.get("formed", None) is not None:
                new_change_dict["broken"] = formed_stereo
            if broken_stereo := change_dict.get("broken", None) is not None:
                new_change_dict["formed"] = broken_stereo
            # raise ValueError(change_dict ,new_change_dict)
            rev_reac.set_atom_stereo_change(atom, **new_change_dict)

        for bond, change_dict in rev_reac._bond_stereo_change.items():
            new_change_dict = {
                stereo_change.value: bond_stereo
                for stereo_change, bond_stereo in change_dict.items()
            }
            if formed_stereo := change_dict.get("formed", None) is not None:
                new_change_dict["broken"] = formed_stereo
            if broken_stereo := change_dict.get("broken", None) is not None:
                new_change_dict["formed"] = broken_stereo
            # raise ValueError(change_dict ,new_change_dict)
            rev_reac.set_bond_stereo_change(bond, **new_change_dict)

        return rev_reac

    def enantiomer(self) -> Self:
        """
        Creates the enantiomer of the StereoCondensedReactionGraph by inversion
        of all chiral stereochemistries. The result can be identical to the
        molecule itself if the molecule is not chiral.

        :return: Enantiomer
        """
        enantiomer = super().enantiomer()
        for atom in self.atoms:
            stereo_change = self.get_atom_stereo_change(atom=atom)
            if stereo_change is not None:
                stereo_change_inverted = {
                    change.value: stereo.invert() if stereo else None
                    for change, stereo in stereo_change.items()
                }
                enantiomer.set_atom_stereo_change(
                    atom=atom, **stereo_change_inverted
                )
        return enantiomer

    def _to_rdmol(
        self, generate_bond_orders=False, charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        
        ts_smg = StereoMolGraph(self) # bond change is now just a bond
        
        for atom, stereo_change_dict in self.atom_stereo_changes.items():
            atom_stereo = next((stereo for stereo_change in (StereoChange.FLEETING,
                                                             StereoChange.BROKEN,
                                                             StereoChange.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if atom_stereo:
                ts_smg.set_atom_stereo(atom, atom_stereo)

        for bond, stereo_change_dict in self.bond_stereo_changes.items():
            bond_stereo = next((stereo for stereo_change in (StereoChange.FLEETING,
                                                             StereoChange.BROKEN,
                                                             StereoChange.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if bond_stereo:
                ts_smg.set_bond_stereo(bond, bond_stereo)

        mol, idx_map_num_dict = ts_smg._to_rdmol(
            generate_bond_orders=False, charge=charge)
        if generate_bond_orders:
            mol = _set_crg_bond_orders(graph=self,
                                 mol=mol,
                                 charge=charge,
                                 idx_map_num_dict=idx_map_num_dict)
        return mol, idx_map_num_dict

    @classmethod
    def from_composed_molgraphs(cls, mol_graphs: Iterable[Self]) -> Self:
        """Creates a MolGraph object from a list of MolGraph objects

        :param mol_graphs: list of MolGraph objects
        :return: Returns Combined MolGraph
        """
        graph = cls(super().from_composed_molgraphs(mol_graphs))
        for mol_graph in mol_graphs:
            graph._atom_stereo_change.update(mol_graph._atom_stereo_change)
            graph._bond_stereo_change.update(mol_graph._bond_stereo_change)
        return graph

    @classmethod
    def from_reactant_and_product_graph(
        cls: type[Self],
        reactant_graph: StereoMolGraph,
        product_graph: StereoMolGraph,
    ) -> Self:
        """Creates a StereoCondensedReactionGraph from reactant and product
        StereoMolGraphs.

        StereoCondensedReactionGraph  is constructed from bond changes from
        reactant to the product. The atoms order and atom types of the reactant
        and product have to be the same.

        :param reactant_graph: reactant of the reaction
        :param product_graph: product of the reaction
        :return: StereoCondensedReactionGraph
        """
        crg = super().from_reactant_and_product_graph(
            reactant_graph, product_graph
        )
        scrg = StereoCondensedReactionGraph(crg)

        all_stereo_atoms = set(reactant_graph._atom_stereo) | set(
            product_graph._atom_stereo
        )

        for atom in all_stereo_atoms:
            r_stereo = reactant_graph.get_atom_stereo(atom)
            p_stereo = product_graph.get_atom_stereo(atom)

            if (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo == p_stereo
            ):
                scrg.set_atom_stereo(atom, r_stereo)

            elif r_stereo is None and p_stereo is not None:
                scrg.set_atom_stereo_change(atom, formed=p_stereo)

            elif p_stereo is None and r_stereo is not None:
                scrg.set_atom_stereo_change(atom, broken=r_stereo)

            elif (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo != p_stereo
            ):
                scrg.set_atom_stereo_change(
                    atom, formed=p_stereo, broken=r_stereo
                )

        all_stereo_bonds = set(reactant_graph._bond_stereo) | set(
            product_graph._bond_stereo
        )

        for bond in all_stereo_bonds:
            r_stereo = reactant_graph.get_bond_stereo(bond)
            p_stereo = product_graph.get_bond_stereo(bond)

            if (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo == p_stereo
            ):
                scrg.set_bond_stereo(bond, r_stereo)

            elif r_stereo is None and p_stereo is not None:
                scrg.set_bond_stereo_change(bond, formed=p_stereo)

            elif p_stereo is None and r_stereo is not None:
                scrg.set_bond_stereo_change(bond, broken=r_stereo)

            elif (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo != p_stereo
            ):
                scrg.set_bond_stereo_change(
                    bond, formed=p_stereo, broken=r_stereo
                )

        for atom in scrg.atoms:
            if (
                (
                    scrg.get_atom_stereo_change(atom) is None
                    or all(
                        stereo is None
                        for stereo in scrg.get_atom_stereo_change(
                            atom
                        ).values()
                    )
                )
                and reactant_graph.get_atom_stereo(atom) is None
                and product_graph.get_atom_stereo(atom) is None
                and len(scrg.bonded_to(atom)) == 4
            ):
                scrg.set_atom_stereo_change(
                    atom,
                    broken=Tetrahedral(scrg.bonded_to(atom), None),
                    formed=Tetrahedral(scrg.bonded_to(atom), None),
                )
                # TODO: add someting for 5 or more substituents

        return scrg

    @classmethod
    def from_reactant_and_product_geometry(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        switching_function: Callable = BondsFromDistance(),
    ) -> Self:
        """Creates a StereoCondensedReactionGraph from reactant and product
        Geometries.

        StereoCondensedReactionGraph  is constructed from bond changes from
        reactant to the product. The atoms order and atom types of the reactant
        and product have to be the same. The switching function is used to
        determine the connectivity of the atoms.

        :param reactant_geo: geometry of the reactant
        :param product_geo: geometry of the product
        :param switching_function: function to define the connectivity from
                                   geometry,
                                   defaults to StepSwitchingFunction()
        :return: StereoCondensedReactionGraph
        """
        reactant_graph = StereoMolGraph.from_geometry(
            reactant_geo, switching_function
        )
        product_graph = StereoMolGraph.from_geometry(
            product_geo, switching_function
        )

        return cls.from_reactant_and_product_graph(
            reactant_graph=reactant_graph, product_graph=product_graph
        )

    @classmethod
    def from_reactant_product_and_ts_geometry(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        ts_geo: Geometry,
        switching_function: Callable = BondsFromDistance(),
    ) -> Self:
        """Creates a StereoCondensedReactionGraph from reactant, product and
        transition state Geometries.

        StereoCondensedReactionGraph  is constructed from bond changes from
        reactant to the product. The atoms order and atom types of the reactant
        and product have to be the same. The switching function is used to
        determine the connectivity of the atoms. Only the stereo information
        is taken from the transition state geometry.

        :param reactant_geo: geometry of the reactant
        :param product_geo: geometry of the product
        :param ts_geo: geometry of the transition state
        :param switching_function: function to define the connectivity from
                                   geometry,
                                   defaults to StepSwitchingFunction()
        :return: CondensedReactionGraph
        """

        crg = CondensedReactionGraph.from_reactant_and_product_geometry(
            reactant_geo=reactant_geo,
            product_geo=product_geo,
            switching_function=switching_function,
        )

        ts_atom_stereo_graph = StereoMolGraph(crg)
        ts_atom_stereo_graph._set_atom_stereo_from_geometry(ts_geo)
        reactant_atom_stereo_graph = StereoMolGraph.from_geometry(
            geo=reactant_geo, switching_function=switching_function
        )
        product_atom_stereo_graph = StereoMolGraph.from_geometry(
            geo=product_geo, switching_function=switching_function
        )

        scrg = cls.from_reactant_and_product_graph(
            reactant_graph=reactant_atom_stereo_graph,
            product_graph=product_atom_stereo_graph,
        )

        for atom in ts_atom_stereo_graph._atom_stereo:
            r_atom_stereo = reactant_atom_stereo_graph.get_atom_stereo(atom)
            ts_atom_stereo = ts_atom_stereo_graph.get_atom_stereo(atom)
            p_atom_stereo = product_atom_stereo_graph.get_atom_stereo(atom)

    
            if (r_atom_stereo is not None
                and p_atom_stereo is not None
                and ts_atom_stereo is not None
                and r_atom_stereo != p_atom_stereo
                and r_atom_stereo != ts_atom_stereo
                and p_atom_stereo != ts_atom_stereo

            ):
                scrg.set_atom_stereo_change(
                    atom=atom,
                    broken=r_atom_stereo,
                    formed=p_atom_stereo,
                    fleeting=ts_atom_stereo,
                )

        return scrg

    def get_isomorphic_mappings(
        self, other: Self, stereo=True, stereo_change=True
    ) -> Iterator[dict[int, int]]:
        """Isomorphic mappings between "self" and "other".

        Generates all isomorphic mappings between "other" and "self".
        All atoms and bonds have to be present in both graphs.
        The Stereochemistry is preserved in the mappings.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """
        return vf2pp_all_isomorphisms(
            self,
            other,
            labels=None,
            color_refine=False, # TODO: implement color refinement
            stereo=stereo,
            stereo_change=stereo_change,
            subgraph=False,
        )

    def get_subgraph_isomorphic_mappings(
        self, other: Self, stereo=True, stereo_change=True
    ) -> Iterator[dict[int, int]]:
        """Subgraph isomorphic mappings from "other" onto "self".

        Generates all node-iduced subgraph isomorphic mappings.
        All atoms of "other" have to be present in "self".
        The bonds of "other" have to be the subset of the bonds of "self"
        relating to the nodes of "other".
        The Stereochemistry is preserved in the mappings.

        :param other: Other Graph to compare with
        :return: Mappings from the atoms of self onto the atoms of other
        :raises TypeError: Not defined for objects different types
        """
        return vf2pp_all_isomorphisms(
            self,
            other,
            labels=None,
            color_refine=False,
            stereo=stereo,
            stereo_change=stereo_change,
            subgraph=True,
        )

    def apply_reaction(
        self,
        reactant: MolGraph,
        mapping: Iterable[dict[int, int]],
        stereo: bool = True
    ) -> Self:
        """
        Applies a reaction to the graph and returns the resulting graph.
        The reactants of the CRG have to be a subgraph of the reactant.
        Mappings from crg atoms to the reactant can be provided.
        Atom numberig from the reactant is kept in the resulting graph.

        :param reaction: Reaction to apply
        :param mapping: Mappings from reactant atoms to crg atoms
        :return: Resulting graph
        """

        scrg = super().apply_reaction(reactant, mapping)

        if stereo:
            for atom, stereo_change_dict in self._atom_stereo_change.items():
                new_change_dict = {}
                for stereo_change, atom_stereo in stereo_change_dict.items():
                    new_atom_stereo = atom_stereo.__class__(
                        atoms = tuple([mapping[a] for a in atom_stereo.atoms]),
                        parity = atom_stereo.parity,)
                    new_change_dict[stereo_change.value] = new_atom_stereo
                scrg.set_atom_stereo_change(mapping[atom], **new_change_dict)

            for bond, stereo_change_dict in self._bond_stereo_change.items():
                new_change_dict = {}
                for stereo_change, bond_stereo in stereo_change_dict.items():
                    new_bond_stereo = bond_stereo.__class__(
                        atoms = tuple([mapping[a] for a in bond_stereo.atoms]),
                        parity = bond_stereo.parity,
                        )
                    new_change_dict[stereo_change.value] = new_bond_stereo
                scrg.set_bond_stereo_change(
                    Bond((mapping[bond[0]], mapping[bond[1]])),
                    **new_change_dict,
                    )

            stereo_change_atoms = [atom for atom in scrg._atom_stereo_change
                                   if atom in scrg._atom_stereo]
            for atom in stereo_change_atoms:
                scrg.delete_atom_stereo(atom)

            stereo_change_bonds = [bond for bond in scrg._bond_stereo_change
                                   if bond in scrg._bond_stereo]
            for bond in stereo_change_bonds:
                scrg.delete_bond_stereo(bond)
                
        return scrg
    

