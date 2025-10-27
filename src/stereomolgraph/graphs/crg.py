from __future__ import annotations

from collections import Counter
from enum import Enum
from typing import TYPE_CHECKING

from stereomolgraph.algorithms.color_refine import (
    numpy_int_multiset_hash,
    label_hash,
    color_refine_crg,
)
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.coords import BondsFromDistance
from stereomolgraph.graph2rdmol import mol_graph_to_rdmol, set_crg_bond_orders
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph

if TYPE_CHECKING:
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
        if self.n_atoms == 0:
            return hash(self.__class__)
        color_array = color_refine_crg(self)
        return int(numpy_int_multiset_hash(color_array))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented

        o_labels = label_hash(other, atom_labels=("atom_type", "reaction"))
        s_labels = label_hash(self, atom_labels=("atom_type", "reaction"))
        o_color_array = color_refine_crg(other, atom_labels=o_labels)
        s_color_array = color_refine_crg(self, atom_labels=s_labels)

        o_colors = {a: int(c) for a,c in zip(other.atoms, o_color_array)}
        s_colors = {a: int(c) for a,c in zip(self.atoms, s_color_array)}

        return any(
                vf2pp_all_isomorphisms(
                    self,
                    other,
                    atom_labels=(s_colors, o_colors),
                    stereo=False,
                    stereo_change=False,
                    subgraph=False,
                )
            )

    def add_bond(self, atom1: int, atom2: int, **attr: Any):
        """
        Adds a bond between atom1 and atom2.

        :param atom1: id of atom1
        :param atom2:   id of atom2
        """
        if "reaction" in attr and not isinstance(attr.get("reaction"), Change):
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
        
    def add_fleeting_bond(self, atom1: int, atom2: int, **attr: Any):
        """
        Adds a bond between atom1 and atom2 with reaction attribute
        set to FLEETING.

        :param atom1: Atom1
        :param atom2: Atom2
        """
        if atom1 in self._atom_attrs and atom2 in self._atom_attrs:
            self.add_bond(atom1, atom2, reaction=Change.FLEETING, **attr)
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
            if (
                self.get_bond_attribute(atom1, atom2, "reaction")
                == Change.FORMED
            ):
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
            if (
                self.get_bond_attribute(atom1, atom2, "reaction")
                == Change.BROKEN
            ):
                b_bonds.add(bond)
        return b_bonds
    
    def get_fleeting_bonds(self) -> set[Bond]:
        """
        Returns all bonds that are fleeting during the reaction.

        :return: fleeting bonds
        """
        f_bonds: set[Bond] = set()
        for bond in self.bonds:
            atom1, atom2 = bond
            if (
                self.get_bond_attribute(atom1, atom2, "reaction")
                == Change.FLEETING
            ):
                f_bonds.add(bond)
        return f_bonds

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
        generate_bond_orders: bool = False,
        allow_charged_fragments: bool = False,
        charge: int = 0,
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        mol, idx_map_num_dict = mol_graph_to_rdmol(
            graph=self,
            generate_bond_orders=False,
            allow_charged_fragments=allow_charged_fragments,
            charge=0,
        )
        set_crg_bond_orders(
            graph=self, mol=mol, idx_map_num_dict=idx_map_num_dict
        )
        return mol, idx_map_num_dict

    def to_rdmol(
        self,
        generate_bond_orders: bool = False,
        allow_charged_fragments: bool = False,
        charge: int = 0,
    ) -> Chem.rdchem.Mol:
        raise NotImplementedError(
            "Rdkit is not able to represent "
            "reactions as condensed reaction graphs."
        )

    def reactant(self, keep_attributes: bool = True) -> MolGraph:
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
                attrs = {"atom_type": self._atom_attrs[atom]["atom_type"]}
            product.add_atom(atom, **attrs)
        for bond in self.bonds:
            bond_reaction = self._bond_attrs[bond].get("reaction", None)
            if bond_reaction is None or bond_reaction == Change.BROKEN:
                if keep_attributes is True:
                    attrs = self._bond_attrs[bond].copy()
                    attrs.pop("reaction", None)
                else:
                    attrs = {}
                product.add_bond(*bond, **attrs)
        return product

    def product(self, keep_attributes: bool = True) -> MolGraph:
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
                attrs = {"atom_type": self._atom_attrs[atom]["atom_type"]}
            product.add_atom(atom, **attrs)
        for bond in self.bonds:
            bond_reaction = self._bond_attrs[bond].get("reaction", None)
            if bond_reaction is None or bond_reaction == Change.FORMED:
                if keep_attributes is True:
                    attrs = self._bond_attrs[bond].copy()
                    attrs.pop("reaction", None)
                else:
                    attrs = {}
                product.add_bond(*bond, **attrs)
        return product

    def _ts(self, keep_attributes: bool = True) -> MolGraph:
        return MolGraph(self)

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
    def from_graphs(
        cls,
        reactant_graph: MolGraph,
        product_graph: MolGraph,
        ts_graph: MolGraph | None = None,
    ) -> Self:
        """Creates a CondensedReactionGraph from reactant and product MolGraphs

        CondensedReactionGraph  is constructed from bond changes from reactant
        to the product. The atoms order and atom types of the reactant and
        product have to be the same.

        :param reactant_graph: reactant of the reaction
        :param product_graph: product of the reaction
        :return: CondensedReactionGraph
        """
        
        r_id_type = set(zip(reactant_graph.atoms, reactant_graph.atom_types))
        p_id_type = set(zip(product_graph.atoms, product_graph.atom_types))
        assert r_id_type == p_id_type
        
        crg = cls()

        for atom, atom_type in zip(reactant_graph.atoms,
                                   reactant_graph.atom_types):
            crg.add_atom(atom,atom_type)

        bonds = set(reactant_graph.bonds) | set(product_graph.bonds)
            
        for bond in bonds:
            if (reactant_graph.has_bond(*bond)
                and product_graph.has_bond(*bond)):
                crg.add_bond(*bond)
            elif reactant_graph.has_bond(*bond):
                crg.add_broken_bond(*bond)
            elif product_graph.has_bond(*bond):
                crg.add_formed_bond(*bond)

        if ts_graph is not None:
            ts_id_type = set(zip(ts_graph.atoms, ts_graph.atom_types))
            assert r_id_type == ts_id_type
            assert p_id_type == ts_id_type

            assert set(ts_graph.bonds).issuperset(bonds), ("TS graph has to "
                        "contain all bonds from reactant and product")

            for bond in ts_graph.bonds:
                if bond not in crg.bonds:
                    crg.add_fleeting_bond(*bond)

        return crg

    @classmethod
    def from_geometries(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        ts_geo: Geometry | None = None,
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
        return cls.from_graphs(reactant_graph=reactant, product_graph=product)
