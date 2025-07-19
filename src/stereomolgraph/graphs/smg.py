from __future__ import annotations

from collections import Counter
from copy import deepcopy
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal, TypeVar

import numpy as np

from stereomolgraph.algorithms.color_refine import color_refine_mg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.coords import are_planar
from stereomolgraph.graph2rdmol import stereo_mol_graph_to_rdmol
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
from stereomolgraph.rdmol2graph import stereo_mol_graph_from_rdmol
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    BondStereo,
    PlanarBond,
    Tetrahedral,
    TrigonalBipyramidal,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping
    from typing import Self

    from rdkit import Chem

    from stereomolgraph.coords import Geometry
    A = TypeVar("A", bound=tuple[int, ...], covariant=True)
    P = TypeVar("P", bound=None | Literal[1, 0, -1], covariant=True)


class StereoMolGraph(MolGraph):
    """
    :class:`MolGraph` with the ability to store stereochemistry information
    for atoms and bonds.

    Two graphs compare equal, if they are isomorphic and have the same
    stereochemistry.
    """
    __slots__ = ("_atom_stereo", "_bond_stereo")

    _atom_stereo: dict[int, AtomStereo]
    _bond_stereo: dict[Bond, BondStereo]

    def __init__(self, mol_graph: None | MolGraph = None):
        super().__init__(mol_graph)
        if mol_graph and isinstance(mol_graph, StereoMolGraph):
            self._atom_stereo = deepcopy(mol_graph._atom_stereo)
            self._bond_stereo = deepcopy(mol_graph._bond_stereo)
        else:
            self._atom_stereo = {}
            self._bond_stereo = {}

    def __hash__(self) -> int:
        color_dict: dict[AtomId, int] = color_refine_mg(self, )
        connectivity_hash = set(Counter(color_dict.values()).items())
        atom_stereo = {s.__class__( tuple([color_dict[a] for a in s.atoms]),
                                   s.parity)
                                   for s in self.atom_stereo.values()}
        bond_stereo = {s.__class__(tuple([color_dict[a] for a in s.atoms]),
                                   s.parity)
                                   for s in self.bond_stereo.values()}
        return hash((frozenset(connectivity_hash),
                     frozenset(atom_stereo),
                     frozenset(bond_stereo)))

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
    ) -> None | AtomStereo:
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

    def set_atom_stereo(self, atom_stereo: AtomStereo):
        """Adds stereo information to the graph

        :param atom: Atoms to be used for chiral information
        :param stereo: Chiral information
        """
        atom = atom_stereo.central_atom
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
    ) -> None | BondStereo:
        """Gets the stereo information of the bond or None
        if it does not exist.
        Raises a ValueError if the bond s not in the graph.

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
        self, bond_stereo: BondStereo
    ):
        """Stets the stereo information of the bond

        :param bond: Bond
        :param bond_stereo: Stereo information of the bond
        """

        bond = Bond(bond_stereo.bond)
        if bond in self._bond_attrs:
            self._bond_stereo[bond] = bond_stereo
        else:
            raise ValueError(f"Bond {bond} is not in the graph")

    def delete_bond_stereo(self, bond: Iterable[int]):
        """Deletes the stereo information of the bond

        :param bond: Bond
        """
        del self._bond_stereo[Bond(bond)]

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
                new_graph.set_atom_stereo(atoms_atom_stereo)

        for _bond, bond_stereo in self._bond_stereo.items():
            if all(atom in atoms for atom in bond_stereo.atoms):
                new_graph.set_bond_stereo(bond_stereo)
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
                enantiomer.set_atom_stereo(stereo.invert())
        return enantiomer

    def _to_rdmol(
        self,
        generate_bond_orders:bool=False,
        allow_charged_fragments:bool=False,
        charge:int=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        """
        Creates a RDKit mol object using the connectivity of the mol graph.
        Stereochemistry is added to the mol object.

        :return: RDKit molecule
        """
        return stereo_mol_graph_to_rdmol(self,
                                          generate_bond_orders=generate_bond_orders,
                                          allow_charged_fragments=allow_charged_fragments,
                                          charge=charge)
        

    @classmethod
    def from_rdmol(cls, rdmol:Chem.Mol, use_atom_map_number:bool=False
                   ) -> Self:
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
        smg = stereo_mol_graph_from_rdmol(cls, rdmol,
                                          use_atom_map_number=use_atom_map_number)
        assert isinstance(smg, cls), ("StereoMolGraph.from_rdmol did not"
                                      " return a StereoMolGraph")
        return smg

    def _set_atom_stereo_from_geometry(self, geo: Geometry):
        for atom in range(geo.n_atoms):
            first_neighbors = list(self.bonded_to(atom))

            # extends the first layer of neighbors to the second layer
            # (if planar)
            # needed to find double bonds
            if len(first_neighbors) < 3:
                pass
            elif (len(first_neighbors) == 3
                  and are_planar(geo.coords[[atom, *first_neighbors]])):
                for neighbor in first_neighbors:
                    next_layer = set(self.bonded_to(neighbor))

                    for i in first_neighbors:
                        next_layer.add(i)
                    next_layer -= set((neighbor, atom))

                    if len(next_layer) != 4:
                        continue

                    elif are_planar(geo.coords.take(tuple(next_layer),
                                                    axis=0)):
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
                        assert len(atom_ids) == 6
                        double_bond = PlanarBond.from_coords(
                            atom_ids, geo.coords.take(atom_ids, axis=0))
                        
                        self.set_bond_stereo(double_bond)

            elif len(first_neighbors) == 3:
                pass
            else:
                first_neighbors_coords = geo.coords[first_neighbors]
                if are_planar(first_neighbors_coords):
                    pass

                elif len(first_neighbors) == 4:
                    stereo_atoms = (atom, *first_neighbors)
                    assert len(stereo_atoms) == 5
                    stereo_coords = geo.coords.take(stereo_atoms, axis=0)
                    atoms_atom_stereo = Tetrahedral.from_coords(stereo_atoms,
                                                                stereo_coords)

                    self.set_atom_stereo(atoms_atom_stereo)

                elif len(first_neighbors) == 5:
                    stereo_atoms = (atom, *first_neighbors)
                    assert len(stereo_atoms) == 6
                    stereo_coords = geo.coords.take(stereo_atoms, axis=0)
                    atoms_atom_stereo = TrigonalBipyramidal.from_coords(
                        stereo_atoms, stereo_coords)
                    self.set_atom_stereo(atoms_atom_stereo)


    @classmethod
    def compose(cls, mol_graphs: Iterable[MolGraph]) -> Self:
        """Creates a MolGraph object from a list of MolGraph objects.
        
        Duplicate nodes or edges are overwritten, such that the resulting
        graph only contains one node or edge with that name. Duplicate
        attributes of duplicate nodes, edges and the stereochemistry are also
        overwritten in order of iteration.

        :param mol_graphs: list of MolGraph objects
        :return: Returns MolGraph
        """

        graph = cls(super().compose(mol_graphs))
        for mol_graph in mol_graphs:
            graph._atom_stereo.update(cls(mol_graph)._atom_stereo)
            graph._bond_stereo.update(cls(mol_graph)._bond_stereo)
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
        self, other: MolGraph, stereo:bool=True
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
            color_refine=True, #TODO: implement color refinement
            stereo=stereo,
            stereo_change=False,
            subgraph=False,
        )
        

    def get_subgraph_isomorphic_mappings(
        self, other: MolGraph, stereo: bool = True
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
        if other.__class__ is not self.__class__:
            return
            yield
        
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

