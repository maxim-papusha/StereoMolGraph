from __future__ import annotations

import warnings
from collections import Counter, defaultdict, deque
from copy import deepcopy
from enum import Enum
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
from stereomolgraph.graphs.mg import AtomId, Bond, MolGraph
from stereomolgraph.graphs.smg import StereoMolGraph
from stereomolgraph.graphs.crg import CondensedReactionGraph
from stereomolgraph import PERIODIC_TABLE, Element
from stereomolgraph.cartesian import are_planar, BondsFromDistance
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.algorithms.color_refine import color_refine_mg
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
        *,
        broken: Optional[AtomStereo] = None,
        fleeting: Optional[AtomStereo] = None,
        formed: Optional[AtomStereo] = None,
    ):
        if 1 != len({(atom := stereo.central_atom) for stereo in
                         (broken, fleeting, formed)
                         if stereo is not None}):
            raise ValueError("Provide stereo information for one atom only")

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
        *,
        broken: Optional[BondStereo] = None,
        fleeting: Optional[BondStereo] = None,
        formed: Optional[BondStereo] = None,

    ):
        if 1 != len({(bond := stereo.bond) for stereo in
                         (broken, fleeting, formed)
                         if stereo is not None}):
            raise ValueError("Provide stereo information for one bond only")
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
                #reactant._bond_stereo[bond] = stereo
                reactant.set_bond_stereo(stereo)

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
                product.set_atom_stereo(stereo)

        for bond, change_dict in self._bond_stereo_change.items():
            if stereo := change_dict[StereoChange.FORMED]:
                product.set_bond_stereo(stereo)

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
            rev_reac.set_atom_stereo_change(**new_change_dict)

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
            rev_reac.set_bond_stereo_change(**new_change_dict)

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
        generate_bond_orders=False,
        allow_charged_fragments=False,
        charge=0
    ) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
        
        ts_smg = StereoMolGraph(self) # bond change is now just a bond
        
        for atom, stereo_change_dict in self.atom_stereo_changes.items():
            atom_stereo = next((stereo for stereo_change in (StereoChange.FLEETING,
                                                             StereoChange.BROKEN,
                                                             StereoChange.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if atom_stereo:
                ts_smg.set_atom_stereo(atom_stereo)

        for bond, stereo_change_dict in self.bond_stereo_changes.items():
            bond_stereo = next((stereo for stereo_change in (StereoChange.FLEETING,
                                                             StereoChange.BROKEN,
                                                             StereoChange.FORMED)
                                if (stereo := stereo_change_dict[stereo_change]) is not None), None)
            if bond_stereo:
                ts_smg.set_bond_stereo(bond_stereo)

        mol, idx_map_num_dict = ts_smg._to_rdmol(
            generate_bond_orders=False,
            allow_charged_fragments=allow_charged_fragments,
            charge=charge)
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
                scrg.set_atom_stereo_change(**new_change_dict)

            for bond, stereo_change_dict in self._bond_stereo_change.items():
                new_change_dict = {}
                for stereo_change, bond_stereo in stereo_change_dict.items():
                    new_bond_stereo = bond_stereo.__class__(
                        atoms = tuple([mapping[a] for a in bond_stereo.atoms]),
                        parity = bond_stereo.parity,
                        )
                    new_change_dict[stereo_change.value] = new_bond_stereo
                scrg.set_bond_stereo_change(**new_change_dict)

            stereo_change_atoms = [atom for atom in scrg._atom_stereo_change
                                   if atom in scrg._atom_stereo]
            for atom in stereo_change_atoms:
                scrg.delete_atom_stereo(atom)

            stereo_change_bonds = [bond for bond in scrg._bond_stereo_change
                                   if bond in scrg._bond_stereo]
            for bond in stereo_change_bonds:
                scrg.delete_bond_stereo(bond)
                
        return scrg
    

