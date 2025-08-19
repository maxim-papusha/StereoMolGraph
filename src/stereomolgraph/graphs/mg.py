from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
from pprint import pformat
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np

from stereomolgraph.algorithms.color_refine import color_refine_mg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.coords import BondsFromDistance
from stereomolgraph.graph2rdmol import mol_graph_to_rdmol
from stereomolgraph.periodic_table import PERIODIC_TABLE, Element
from stereomolgraph.rdmol2graph import mol_graph_from_rdmol
from stereomolgraph.xyz2graph import connectivity_from_geometry

if TYPE_CHECKING:
    from collections.abc import (
        Collection,
        Iterable,
        Mapping,
        Sequence,
    )
    from typing import Any, Optional, Self, TypeAlias, TypeVar

    from rdkit import Chem  # type: ignore

    from stereomolgraph.coords import Geometry

    N = TypeVar("N", bound=int,)

AtomId: TypeAlias = int

Bond: TypeAlias = frozenset[AtomId]


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
    ) -> Collection[AtomId]:
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
    def neighbors(
        self,
    ) -> Mapping[AtomId, set[AtomId]]:
        """
        :return: Returns all neighbors of the atoms in the MolGraph
        """
        return MappingProxyType(self._neighbors)

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
        self, atom: AtomId, atom_type: int | str | Element, **attr: Any
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
                self.remove_bond(atom, n)

    def get_atom_attribute(self, atom: AtomId, attr: str) -> Optional[Any]:
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
        return Bond({atom1, atom2}) in self._bond_attrs

    def add_bond(self, atom1: AtomId, atom2: AtomId, **attr: Any):
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
        self._neighbors[atom1].discard(atom2)
        self._neighbors[atom2].discard(atom1)

    def get_bond_attribute(
        self,
        atom1: AtomId,
        atom2: AtomId,
        attr: str,
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

    def delete_bond_attribute(self, atom1: AtomId, atom2: AtomId, attr: str):
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

    def connectivity_matrix(
        self,
    ) -> np.ndarray[tuple[int, int], np.dtype[np.int8]]:
        """
        Returns a connectivity matrix of the graph as a list of lists.
        Order is the same as in self.atoms()
        1 if nodes are connected, 0 if not.

        :return: Connectivity matrix as list of lists
        """
        n = len(self.atoms)
        matrix = np.zeros((n, n), dtype=np.int8)
        atomid_index_dict = {id: index for index, id in enumerate(self.atoms)}

        for a1, a2 in self.bonds:
            matrix[atomid_index_dict[a1]][atomid_index_dict[a2]] = 1
            matrix[atomid_index_dict[a2]][atomid_index_dict[a1]] = 1
        return matrix

    def _to_rdmol(
        self,
        generate_bond_orders: bool = False,
        allow_charged_fragments: bool = False,
        charge: int = 0,
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        return mol_graph_to_rdmol(
            self,
            generate_bond_orders=generate_bond_orders,
            allow_charged_fragments=allow_charged_fragments,
            charge=charge,
        )

    def to_rdmol(
        self,
        generate_bond_orders: bool = False,
        allow_charged_fragments: bool = False,
        charge: int = 0,
    ) -> Chem.rdchem.Mol:
        mol, _ = self._to_rdmol(
            generate_bond_orders=generate_bond_orders,
            allow_charged_fragments=allow_charged_fragments,
            charge=charge,
        )
        for atom in mol.GetAtoms():  # type: ignore
            atom.SetAtomMapNum(0, strict=True)  # type: ignore
        return mol

    @classmethod
    def from_rdmol(
        cls, rdmol: Chem.Mol, use_atom_map_number: bool = False
    ) -> Self:
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
        return mol_graph_from_rdmol(
            cls, rdmol, use_atom_map_number=use_atom_map_number
        )  # type: ignore

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

    def node_connected_component(self, atom: int) -> set[AtomId]:
        """
        :param atom: atom id
        :return: Returns the connected component that includes atom_id
        """
        visited: set[AtomId] = set()
        stack = [atom]
        while stack:
            node = stack.pop()

            if node not in visited:
                visited.add(node)
            for neighbor in self.bonded_to(node):
                if neighbor not in visited:
                    stack.append(neighbor)
        return visited

    def connected_components(self) -> list[set[int]]:
        """
        :return: Returns the connected components of the graph
        """
        visited: set[AtomId] = set()
        components: list[set[int]] = []

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
        matrix: np.ndarray,
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
    def compose(cls, mol_graphs: Iterable[MolGraph]) -> Self:
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
        threshold: float = 0.5,
        include_bond_order: bool = False,
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
    def from_geometry_and_bond_order_matrix(
        cls,
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
        cls,
        geo: Geometry,
        switching_function: BondsFromDistance = BondsFromDistance(),
    ) -> Self:
        return connectivity_from_geometry(cls, geo, switching_function)

    def color_refine(self) -> Mapping[AtomId, int]:
        """
        Color refines the graph by assigning a color to each atom based on its
        type and connectivity. The color is a unique integer for each unique
        atom type and connectivity pattern.

        :return: A dictionary mapping atom ids to their color labels.
        """
        return color_refine_mg(self)

    def __hash__(self) -> int:
        return hash(frozenset(Counter(self.color_refine().values()).items()))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            other_atom_labels = other.color_refine()
            self_atom_labels = self.color_refine()

            return any(
                vf2pp_all_isomorphisms(
                    self,
                    other,
                    atom_labels=(self_atom_labels, other_atom_labels),
                    stereo=False,
                    stereo_change=False,
                    subgraph=False,
                )
            )
        else:
            return NotImplemented

    def is_isomorphic(self, other: Self) -> bool:
        return self == other

    def __str__(self) -> str:
        a_list = sorted(
            (a, a_type.symbol)
            for a, a_type in zip(self.atoms, self.atom_types)
        )
        b_list = sorted(tuple(sorted(bond)) for bond in self.bonds)

        pretty_str = pformat(
            [["Atoms", a_list], ["Bonds", b_list]],
            indent=0,
            width=120,
            compact=True,
            sort_dicts=True,
        )
        return f"{self.__class__.__name__}\n{pretty_str}".translate(
            str.maketrans("", "", ",\"'[]")
        )

    def _ipython_display_(self) -> None:
        print(self.__repr__())
