from __future__ import annotations


from enum import Enum

from typing import TYPE_CHECKING

import numpy as np
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
from stereomolgraph.cartesian import BondsFromDistance
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms


from stereomolgraph.graph2rdmol import (
    _mol_graph_to_rdmol,
    _set_crg_bond_orders
)

if TYPE_CHECKING:

    from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
    from typing import Any, Optional, TypeAlias, TypeVar
    
    from rdkit import Chem

    from stereomolgraph.cartesian import Geometry

    G = TypeVar("G", bound="CondensedReactionGraph", covariant=True)



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
            

    def get_formed_bonds(self) -> set[Bond]:
        """
        Returns all bonds that are formed during the reaction

        :return: formed bonds
        """
        return {
            bond
            for bond in self.bonds
            if self.get_bond_attribute(*bond, "reaction") == BondChange.FORMED
        }

    def get_broken_bonds(self) -> set[Bond]:
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
        self,
        generate_bond_orders=False,
        allow_charged_fragments=False,
        charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        mol, idx_map_num_dict = _mol_graph_to_rdmol(graph=self,
                                                    generate_bond_orders=False,
                                                    allow_charged_fragments=allow_charged_fragments,
                                                    charge=0)
        _set_crg_bond_orders(graph=self, mol=mol, idx_map_num_dict=idx_map_num_dict)
        return mol, idx_map_num_dict

    def to_rdmol(self, *args) -> Chem.rdchem.RWMol:
        raise NotImplementedError(
            "Rdkit is not able to represent "
            "reactions as condensed reaction graphs."
        )

    def reactant(self, keep_attributes:bool=True) -> MolGraph:
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

    def product(self, keep_attributes:bool=True) -> MolGraph:
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

    def reverse_reaction(self) -> G:
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
        cls: type[G], reactant_graph: MolGraph, product_graph: MolGraph
    ) -> G:
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
    ) -> G:
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

    def get_isomorphic_mappings(self, other: G) -> Iterator[dict[int, int]]:
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
    ) -> G:
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

