from __future__ import annotations

from collections import Counter
from enum import Enum
from typing import TYPE_CHECKING

from stereomolgraph.algorithms.color_refine import color_refine_mg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.coords import BondsFromDistance
from stereomolgraph.graph2rdmol import mol_graph_to_rdmol, set_crg_bond_orders
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph

if TYPE_CHECKING:
    
    # Self is included in typing from 3.11

    from collections.abc import Iterator
    from typing import Any, Self

    from rdkit import Chem

    from stereomolgraph.coords import Geometry


class Change(Enum):
    FORMED = "formed"
    FLEETING = "fleeting"
    BROKEN = "broken"

    def __repr__(self) -> str:
        return self.name


class CondensedReactionGraph(MolGraph):
    """
    Graph representing a reaction. Atoms are nodes and (potentially changing)
    bonds are edges. Every node has to have an attribute "atom_type" of type
    Element. Edges can have an attribute "reaction" of type Change.
    This is used to represent the change in connectivity during the reaction.

    Two graphs are equal, iff. they are isomporhic and of the same type.
    """
    __slots__: tuple[str, ...] = tuple()
    _atom_attrs: dict[AtomId, dict[str, Any]]
    _neighbors: dict[AtomId, set[AtomId]]
    _bond_attrs: dict[Bond, dict[str, Any]]

    def __hash__(self) -> int:
        r_color_dict = color_refine_mg(self.reactant())
        p_color_dict = color_refine_mg(self.product())
        ts_colors = color_refine_mg(self)
        color_dict = {a: (r_color_dict[a], ts_colors[a], p_color_dict[a])
                      for a in self.atoms}
        return hash(frozenset(Counter(color_dict.values()).items()))

    def add_bond(self, atom1: int, atom2: int, **attr:Any):
        """
        Adds a bond between atom1 and atom2.

        :param atom1: id of atom1
        :param atom2:   id of atom2
        """
        if "reaction" in attr and not isinstance(
            attr.get("reaction"), Change
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
        if attr == "reaction" and not isinstance(value, Change):
            raise ValueError("reaction bond has to have reaction attribute")
        super().set_bond_attribute(atom1, atom2, attr, value)

    def add_formed_bond(self, atom1: int, atom2: int, **attr: Any):
        """
        Adds a bond between atom1 and atom2 with reaction attribute
        set to FORMED.

        :param atom1: Atom1
        :param atom2: Atom2
        """

        attr["reaction"] = Change.FORMED
        if atom1 in self._atom_attrs and atom2 in self._atom_attrs:
            self.add_bond(atom1, atom2, **attr)
        else:
            raise ValueError("Atoms have to be in the graph")

    def add_broken_bond(self, atom1: int, atom2: int, **attr: Any):
        """
        Adds a bond between atom1 and atom2 with reaction attribute
        set to BROKEN.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        if atom1 in self._atom_attrs and atom2 in self._atom_attrs:
            self.add_bond(atom1, atom2, reaction=Change.BROKEN, **attr)
        else:
            raise ValueError("Atoms have to be in the graph")
            

    def get_formed_bonds(self) -> set[Bond]:
        """
        Returns all bonds that are formed during the reaction

        :return: formed bonds
        """
        f_bonds: set[Bond] = set()
        for bond in self.bonds:
            atom1, atom2 = bond
            if self.get_bond_attribute(atom1, atom2, "reaction") == Change.FORMED:
                f_bonds.add(bond)
        return f_bonds

    def get_broken_bonds(self) -> set[Bond]:
        """
        Returns all bonds that are broken during the reaction

        :return: broken bonds
        """
        b_bonds: set[Bond] = set()
        for bond in self.bonds:
            atom1, atom2 = bond
            if self.get_bond_attribute(atom1, atom2, "reaction") == Change.BROKEN:
                b_bonds.add(bond)
        return b_bonds
        
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

    def _to_rdmol(
        self,
        generate_bond_orders:bool=False,
        allow_charged_fragments:bool=False,
        charge:int=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        mol, idx_map_num_dict = mol_graph_to_rdmol(graph=self,
                                                    generate_bond_orders=False,
                                                    allow_charged_fragments=allow_charged_fragments,
                                                    charge=0)
        set_crg_bond_orders(graph=self, mol=mol, idx_map_num_dict=idx_map_num_dict)
        return mol, idx_map_num_dict

    def to_rdmol(
        self,
        generate_bond_orders:bool=False,
        allow_charged_fragments:bool = False,
        charge:int=0
    ) -> Chem.rdchem.Mol:
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
                bond_reaction is None or bond_reaction == Change.BROKEN
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
                bond_reaction is None or bond_reaction == Change.FORMED
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
            if bond_reaction == Change.FORMED:
                rev_reac.add_broken_bond(*bond)
            elif bond_reaction == Change.BROKEN:
                rev_reac.add_formed_bond(*bond)

        return rev_reac

    @classmethod
    def from_reactant_and_product_graph(
        cls, reactant_graph: MolGraph, product_graph: MolGraph
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
                atom_type=reactant_graph.get_atom_type(atom),
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
                crg.add_bond(*bond, reaction=Change.BROKEN)
            elif product_graph.has_bond(*bond):
                crg.add_bond(*bond, reaction=Change.FORMED)
        return crg

    @classmethod
    def from_reactant_and_product_geometry(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        switching_function: BondsFromDistance = BondsFromDistance(),
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

    def get_isomorphic_mappings(self, other: MolGraph) -> Iterator[dict[int, int]]:
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
            color_refine=True, # TODO: implement color refinement
            stereo=False,
            stereo_change=False,
            subgraph=False,
            # labels=["reaction"],
        )



