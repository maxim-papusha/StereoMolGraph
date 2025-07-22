from __future__ import annotations

import sys

from collections import defaultdict, Counter
from copy import deepcopy
from types import MappingProxyType
from typing import TYPE_CHECKING, Generic

from stereomolgraph.algorithms.color_refine import color_refine_mg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.coords import BondsFromDistance
from stereomolgraph.graph2rdmol import set_crg_bond_orders
from stereomolgraph.graphs.crg import CondensedReactionGraph, Change
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
from stereomolgraph.graphs.smg import StereoMolGraph
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    BondStereo,
    Stereo,
    Tetrahedral,
)

if TYPE_CHECKING:

    import sys
    from collections.abc import Iterable, Iterator, Mapping
    from typing import Optional
    
    from rdkit import Chem

    from stereomolgraph.coords import Geometry
       
# Self is included in typing from 3.11
if sys.version_info >= (3, 11):
    from typing import Self, TypeVar
else:
    from typing_extensions import Self, TypeVar

S = TypeVar("S", bound="Stereo", contravariant=True)


class ChangeDict(dict[Change, None | S], Generic[S]):
    def __missing__(self, key: Change) -> None:
        if key in Change:
            return None
        else:
            raise KeyError(f"{key} not in {self.__class__.__name__}")
    

class StereoCondensedReactionGraph(StereoMolGraph, CondensedReactionGraph):
    """
    :class:`CondenedReactionGraph` with the ability to store stereochemistry
    information for atoms and (potentially changing) bonds.
    """

    __slots__ = ("_atom_stereo_change", "_bond_stereo_change")
    _atom_stereo_change: defaultdict[AtomId, ChangeDict[AtomStereo]]
    _bond_stereo_change: defaultdict[Bond, ChangeDict[BondStereo]]

    def __init__(self, mol_graph: Optional[MolGraph] = None):
        super().__init__(mol_graph)
        self._atom_stereo_change = defaultdict(ChangeDict[AtomStereo])
        self._bond_stereo_change = defaultdict(ChangeDict[BondStereo])
        
        if mol_graph and isinstance(mol_graph, StereoCondensedReactionGraph):
            self._atom_stereo_change.update(mol_graph._atom_stereo_change)
            self._bond_stereo_change.update(mol_graph._bond_stereo_change)
    
    def __hash__(self) -> int:
        # Get molecular graphs and their colorings
        r = self.reactant()
        p = self.product()
        r_colors = color_refine_mg(r)
        p_colors = color_refine_mg(p)
        ts_colors = color_refine_mg(self)
        colors = {a : hash((r_colors[a], p_colors[a], ts_colors[a]))
                  for a in self.atoms}
        
        # Helper function to compute stereo hashes
        def stereo_hash(stereo):
            return hash(stereo.__class__(
                tuple(colors[a] for a in stereo.atoms
            ), stereo.parity))
        
        # Compute atom and bond stereo hashes
        atom_stereo = {a: stereo_hash(s) for a, s in self.atom_stereo.items()}
        bond_stereo = {b: stereo_hash(s) for b, s in self.bond_stereo.items()}
        
        # Compute stereo change hashes
        def stereo_change_hash(change_dict):
            return tuple(
                (change, stereo_hash(s))
                for change, s in change_dict.items()
                if s is not None
            )
        
        atom_stereo_change = {
            a: hash(stereo_change_hash(change_dict))
            for a, change_dict in self._atom_stereo_change.items()
        }
        
        bond_stereo_change = {
            b: hash(stereo_change_hash(change_dict))
            for b, change_dict in self._bond_stereo_change.items()
        }
        
        # Combine all hashes
        components = (
            colors,
            atom_stereo,
            bond_stereo,
            atom_stereo_change,
            bond_stereo_change
        )
        
        return hash(tuple(
            hash(frozenset(Counter(d.values()).items()))
            for d in components
        ))

    @property
    def atom_stereo_changes(self) -> Mapping[AtomId, ChangeDict[AtomStereo]]:
        return MappingProxyType(self._atom_stereo_change)

    @property
    def bond_stereo_changes(self) -> Mapping[Bond, ChangeDict[BondStereo]]:
        return MappingProxyType(self._bond_stereo_change)

    def get_atom_stereo_change(
        self, atom: int
    ) -> None | Mapping[Change, AtomStereo | None]:
        
        if atom in self._atom_attrs:
            if atom in self._atom_stereo_change:
                return MappingProxyType(self._atom_stereo_change[atom])
            else:
                return None
        else:
            raise ValueError(f"Atom {atom} not in graph")

    def get_bond_stereo_change(
        self, bond: Iterable[int]
        ) -> None | Mapping[Change, BondStereo | None]:
        
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
        *,
        broken: Optional[AtomStereo] = None,
        fleeting: Optional[AtomStereo] = None,
        formed: Optional[AtomStereo] = None,
    ):
        atoms: set[int] = set()
        for stereo in (broken, fleeting, formed):
            if stereo is not None:
                atoms.add(stereo.central_atom)
        if len(atoms) != 1:
            raise ValueError("Provide stereo information for one atom only")

        if (atom := atoms.pop()) not in self._atom_attrs:
            raise ValueError(f"Atom {atom} not in graph")
        
        self._atom_stereo_change[atom] = ChangeDict[AtomStereo]()
        for stereo_change, atom_stereo in {
            Change.BROKEN: broken,
            Change.FLEETING: fleeting,
            Change.FORMED: formed,
        }.items():
            if atom_stereo:
                self._atom_stereo_change[atom][stereo_change] = atom_stereo

    def set_bond_stereo_change(
        self,
        *,
        broken: Optional[BondStereo] = None,
        fleeting: Optional[BondStereo] = None,
        formed: Optional[BondStereo] = None,

    ):
        bonds: set[frozenset[int]] = set()
        for stereo in (broken, fleeting, formed):
            if stereo is not None:
                bonds.add(Bond(stereo.bond))
        if len(bonds) != 1:
            raise ValueError("Provide stereo information for one atom only")
        if (bond := bonds.pop()) not in self._bond_attrs:
            raise ValueError(f"Bond {bond} not in graph")
        
        self._bond_stereo_change[bond] = ChangeDict[BondStereo]()
        for stereo_change, bond_stereo in {
            Change.BROKEN: broken,
            Change.FORMED: formed,
            Change.FLEETING: fleeting
        }.items():
            if bond_stereo:
                self._bond_stereo_change[bond][stereo_change] = bond_stereo

    def delete_atom_stereo_change(
        self, atom: AtomId, stereo_change: Optional[Change] = None
    ):
        if stereo_change is None:
            del self._atom_stereo_change[atom]
        else:
            del self._atom_stereo_change[atom][stereo_change]

    def delete_bond_stereo_change(
        self, bond: Iterable[AtomId],
        stereo_change: Optional[Change] = None
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
            
        for _atom, stereo_change in self.atom_stereo_changes.items():
            for _change, stereo in stereo_change.items():
                if stereo is not None:
                    active_atoms.update(stereo.atoms)
        for _bond, stereo_change in self.bond_stereo_changes.items():
            for _change, stereo in stereo_change.items():
                if stereo is not None:
                    active_atoms.update(stereo.atoms)

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

        atom_stereo_change:defaultdict[AtomId, ChangeDict[AtomStereo]] = defaultdict(ChangeDict[AtomStereo])
        
        for atom, stereo_change_dict in self._atom_stereo_change.items():
            for stereo_change, atom_stereo in stereo_change_dict.items():
                if atom_stereo is None:
                    continue
                new_stereo = atom_stereo.__class__(tuple(mapping.get(atom, atom) for atom in atom_stereo.atoms),
                        atom_stereo.parity)
                atom_stereo_change[mapping[atom]][stereo_change] = new_stereo

        bond_stereo_change: defaultdict[Bond, ChangeDict[BondStereo]] = defaultdict(ChangeDict[BondStereo])
        
        for bond, stereo_change_dict in self._bond_stereo_change.items():
            for stereo_change, bond_stereo in stereo_change_dict.items():
                if bond_stereo is None:
                    continue
                new_bond = Bond(mapping[a] for a in bond)
                new_stereo = bond_stereo.__class__(
                    tuple(mapping.get(atom, atom)
                                  for atom in bond_stereo.atoms),
                    bond_stereo.parity,
                )
                bond_stereo_change[new_bond][stereo_change] = new_stereo
            
        if copy is True:
            relabeled_scrg._atom_stereo_change = atom_stereo_change
            relabeled_scrg._bond_stereo_change = bond_stereo_change
        else:
            self._atom_stereo_change = atom_stereo_change
            self._bond_stereo_change = bond_stereo_change

        return relabeled_scrg

    def reactant(self, keep_attributes:bool=True) -> StereoMolGraph:
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
            if stereo := change_dict[Change.BROKEN]:
                reactant._atom_stereo[atom] = stereo
                

        for _bond, change_dict in self._bond_stereo_change.items():
            if stereo := change_dict[Change.BROKEN]:
                #reactant._bond_stereo[bond] = stereo
                reactant.set_bond_stereo(stereo)

        return reactant

    def product(self, keep_attributes:bool=True) -> StereoMolGraph:
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

        for _atom, change_dict in self._atom_stereo_change.items():
            if stereo := change_dict[Change.FORMED]:
                product.set_atom_stereo(stereo)

        for _bond, change_dict in self._bond_stereo_change.items():
            if stereo := change_dict[Change.FORMED]:
                product.set_bond_stereo(stereo)

        return product

    def reverse_reaction(self) -> Self:
        """Creates the reaction in the opposite direction.

        Broken bonds and stereochemistry changes are turned into formed
        and the other way around.

        :return: Reversed reaction
        """
        rev_reac = super().reverse_reaction()
        
        for _atom, atom_change_dict in rev_reac._atom_stereo_change.items():
            new_atom_change_dict = {"fleeting": atom_change_dict[Change.FLEETING],
                               "broken": atom_change_dict[Change.FORMED],
                               "formed": atom_change_dict[Change.BROKEN]}

            rev_reac.set_atom_stereo_change(**new_atom_change_dict)

        for _bond, bond_change_dict in rev_reac._bond_stereo_change.items():
            new_bond_change_dict = {"fleeting": bond_change_dict[Change.FLEETING],
                               "broken": bond_change_dict[Change.FORMED],
                               "formed": bond_change_dict[Change.BROKEN]}
            rev_reac.set_bond_stereo_change(**new_bond_change_dict)
            
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
                    **stereo_change_inverted
                )
        return enantiomer

    def _to_rdmol(
        self,
        generate_bond_orders:bool=False,
        allow_charged_fragments:bool=False,
        charge:int=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        
        ts_smg = StereoMolGraph(self) # bond change is now just a bond
        
        for _atom, stereo_change_dict in self.atom_stereo_changes.items():
            atom_stereo = next((stereo for stereo_change in (Change.FLEETING,
                                                             Change.BROKEN,
                                                             Change.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if atom_stereo:
                ts_smg.set_atom_stereo(atom_stereo)

        for _bond, stereo_change_dict in self.bond_stereo_changes.items():
            bond_stereo = next((stereo for stereo_change in (Change.FLEETING,
                                                             Change.BROKEN,
                                                             Change.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if bond_stereo:
                ts_smg.set_bond_stereo(bond_stereo)

        mol, idx_map_num_dict = ts_smg._to_rdmol(
            generate_bond_orders=False,
            allow_charged_fragments=allow_charged_fragments,
            charge=charge)
        if generate_bond_orders:
            mol = set_crg_bond_orders(graph=self,
                                 mol=mol,
                                 charge=charge,
                                 idx_map_num_dict=idx_map_num_dict)
        return mol, idx_map_num_dict

    @classmethod
    def compose(cls, mol_graphs: Iterable[MolGraph]) -> Self:
        """Creates a MolGraph object from a list of MolGraph objects

        :param mol_graphs: list of MolGraph objects
        :return: Returns Combined MolGraph
        """
        graph = cls(super().compose(mol_graphs))
        for mol_graph in mol_graphs:
            graph._atom_stereo_change.update(cls(mol_graph)._atom_stereo_change)
            graph._bond_stereo_change.update(cls(mol_graph)._bond_stereo_change)
        return graph

    @classmethod
    def from_reactant_and_product_graph(
        cls,
        reactant_graph: MolGraph,
        product_graph: MolGraph,
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
        scrg = cls(crg)

        reactant_graph = StereoMolGraph(reactant_graph)
        product_graph = StereoMolGraph(product_graph)

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
                scrg.set_atom_stereo(r_stereo)

            elif r_stereo is None and p_stereo is not None:
                scrg.set_atom_stereo_change(formed=p_stereo)

            elif p_stereo is None and r_stereo is not None:
                scrg.set_atom_stereo_change(broken=r_stereo)

            elif (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo != p_stereo
            ):
                scrg.set_atom_stereo_change(
                    formed=p_stereo, broken=r_stereo
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
                scrg.set_bond_stereo(r_stereo)

            elif r_stereo is None and p_stereo is not None:
                scrg.set_bond_stereo_change(formed=p_stereo)

            elif p_stereo is None and r_stereo is not None:
                scrg.set_bond_stereo_change(broken=r_stereo)

            elif (
                r_stereo is not None
                and p_stereo is not None
                and r_stereo != p_stereo
            ):
                scrg.set_bond_stereo_change(
                    formed=p_stereo, broken=r_stereo
                )

        for atom in scrg.atoms:
            isnon = (
                    scrg.get_atom_stereo_change(atom) is None
                    or (change_d := scrg.get_atom_stereo_change(atom)) is None
                    or all (stereo is None
                        for stereo in change_d.values()))
            if (isnon
                and reactant_graph.get_atom_stereo(atom) is None
                and product_graph.get_atom_stereo(atom) is None
                and len(scrg.bonded_to(atom)) == 4
            ):
                t_atoms = (atom, *scrg.bonded_to(atom))
                assert len(t_atoms) == 5
                scrg.set_atom_stereo_change(
                    broken=Tetrahedral(t_atoms, None),
                    formed=Tetrahedral(t_atoms, None),
                )
                # TODO: add someting for 5 or more substituents

        return scrg

    @classmethod
    def from_reactant_and_product_geometry(
        cls,
        reactant_geo: Geometry,
        product_geo: Geometry,
        switching_function: BondsFromDistance = BondsFromDistance(),
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
        switching_function: BondsFromDistance = BondsFromDistance(),
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
                    broken=r_atom_stereo,
                    formed=p_atom_stereo,
                    fleeting=ts_atom_stereo,
                )

        return scrg

    def get_isomorphic_mappings(
        self, other: MolGraph, stereo:bool=True, stereo_change:bool=True
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
            color_refine=True, # TODO: implement color refinement
            stereo=stereo,
            stereo_change=stereo_change,
            subgraph=False,
        )

    def get_subgraph_isomorphic_mappings(
        self, other: MolGraph, stereo:bool=True, stereo_change:bool=True
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

    

