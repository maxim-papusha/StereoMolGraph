from __future__ import annotations

import itertools
from collections import deque
from typing import TYPE_CHECKING

from stereomolgraph import StereoMolGraph, StereoCondensedReactionGraph
from stereomolgraph.graph import StereoChange

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable
    from typing import Optional
    from stereomolgraph.graph import (
         AtomId, 
         AtomStereo, 
         BondStereo,
         _BaseAchiralStereo,
         _BaseChiralStereo)


def generate_stereoisomers(graph: StereoMolGraph,
                           enantiomers:bool =  True,
                           atoms: Optional[AtomId] = None,
                           bonds: Optional[Iterable[AtomId, AtomId]] = None) -> set[StereoMolGraph]:
        """Generates all stereoisomers of a StereoMolGraph by generation of
        all combinations of parities. Only includes stereocenter which have a
        parity of None. If a parity is set, it is not changed.
        
        If include_enantiomers is True, both enantiomers of a stereoisomer are
        included, if it is False, only one enantiomer is included.        

        :param enantiomers: If True, both enantiomers are included,
        :param: sets if both enantiomers should be included, default: Ture
        :return: All possible stereoisomers
        """
        if atoms is None:
             atom_stereos = (stereo.get_isomers() for a in graph.atoms
                             if ((stereo := graph.get_atom_stereo(a))
                                 and stereo.parity is None))
        else:
            atom_stereos = (stereo.get_isomers() for a in atoms
                            if (stereo := graph.get_atom_stereo(a) ) is not None)
        
        if bonds is None:
            bond_stereos = (stereo.get_isomers() for b in graph.bonds
                            if ((stereo := graph.get_bond_stereo(b)) and stereo.parity is None))
        else:
            bond_stereos = (stereo.get_isomers() for b in bonds
                            if (stereo := graph.get_bond_stereo(b)) is not None)

        isomers = set()
        enantiomers_set = set()

        for a_stereos, b_stereos in itertools.product(itertools.product(*atom_stereos),
                                                      itertools.product(*bond_stereos)):
            stereoisomer = graph.copy()
            for a_stereo in a_stereos:
                stereoisomer.set_atom_stereo(a_stereo.atoms[0], a_stereo)
            for b_stereo in b_stereos:
                stereoisomer.set_bond_stereo(b_stereo.atoms[2:4], b_stereo)
            
            if stereoisomer not in enantiomers_set:
                isomers.add(stereoisomer)
            
                if not enantiomers:
                    enantiomers_set.add(stereoisomer.enantiomer())

        return isomers
        

def generate_fleeting_stereoisomers(graph: StereoCondensedReactionGraph,
                               enantiomers:bool =  True,
                               atoms: Optional[AtomId] = None,
                               bonds: Optional[Iterable[AtomId, AtomId]] = None
                               ) -> set[StereoCondensedReactionGraph]:
        # TODO: extend to more than fleeting stereochemistry 
        # add checks if fleeting stereochemistry is valid relative to formed and broken

        if atoms is None:
             atom_stereos = [stereo.get_isomers() for a in graph.atoms
                             if ((stereo_change_dict := graph.get_atom_stereo_change(a))
                                 and (stereo := stereo_change_dict[StereoChange.FLEETING]))
                                 and stereo.parity is None]
        else:
            atom_stereos = [stereo.get_isomers() for a in atoms
                             if ((stereo_change_dict := graph.get_atom_stereo_change(a))
                                 and (stereo := stereo_change_dict[StereoChange.FLEETING]))
                                 and stereo.parity is None]
        
        if bonds is None:
            bond_stereos = [stereo.get_isomers() for b in graph.bonds
                             if ((stereo_change_dict := graph.get_bond_stereo_change(b))
                                 and (stereo := stereo_change_dict[StereoChange.FLEETING]))
                                 and stereo.parity is None]
        else:
            bond_stereos = [stereo.get_isomers() for b in bonds
                             if ((stereo_change_dict := graph.get_bond_stereo_change(b))
                                 and (stereo := stereo_change_dict[StereoChange.FLEETING]))
                                 and stereo.parity is None]

        isomers = []
        enantiomers_set = set()

        for a_stereos, b_stereos in itertools.product(itertools.product(*atom_stereos),
                                                      itertools.product(*bond_stereos)):
            stereoisomer = graph.copy()
            for a_stereo in a_stereos:
                stereoisomer.set_atom_stereo_change(a_stereo.atoms[0], fleeting=a_stereo)
            for b_stereo in b_stereos:
                stereoisomer.set_bond_stereo_change(b_stereo.atoms[2:4], fleeting=b_stereo)
            
            if stereoisomer not in enantiomers_set:
                isomers.append(stereoisomer)
            
                if not enantiomers:
                    enantiomers_set.add(stereoisomer.enantiomer())

        return isomers
                

def old_generate_stereoisomers(graph: G, enantiomers=True, ) -> Collection[G]:

        isomers = []

        stereo_nones: list[
            tuple[int | frozenset[int], AtomStereo | BondStereo]
        ] = []

        for atom_bond, atom_bund_stereo in graph.stereo.items():
            if atom_bund_stereo.parity is None:
                stereo_nones.append((atom_bond, atom_bund_stereo))

        combs: list[
            list[tuple[int | frozenset[int], AtomStereo | BondStereo]]
        ] = [[]]

        for atom_bond, stereo in stereo_nones:
            old_combs, combs = combs, []

            if isinstance(stereo, _BaseAchiralStereo):
                parities = (0,)
            elif isinstance(stereo, _BaseChiralStereo):
                parities = (-1, 1)
            for comb, parity in itertools.product(
                old_combs, parities
            ):
                combs.append(
                    [
                        (
                            atom_bond,
                            stereo.__class__(stereo.atoms, parity),
                        ),
                        *comb,
                    ]
                )

        for comb in combs:
            isomer = graph.copy()
            for atom_bond, stereo in comb:
                if isinstance(atom_bond, int):
                    isomer.set_atom_stereo(atom_bond, stereo)
                elif isinstance(atom_bond, frozenset):
                    isomer.set_bond_stereo(atom_bond, stereo)
            if isomer not in isomers:
                if enantiomers is True:
                    isomers.append(isomer)
                elif enantiomers is False:
                    if isomer.enantiomer() not in isomers:
                        isomers.append(isomer)
        return isomers


def topological_symmetry_number(graph: StereoMolGraph) -> int:
        """
        Calculated from the number of graph isomorphisms which conserve the
        stereo information.
        symmetry_number = internal_symmetry_number * rotational_symmetry_number
        TODO: add paper reference
        """

        if any(stereo.parity is None for stereo in graph.stereo.values()):
            raise NotImplementedError(
                "all stereocenters have to be defined"
                " to calculate the symmetry number"
            )

        mappings = graph.get_isomorphic_mappings(graph)
        return deque(enumerate(mappings, 1), maxlen=1)[0][0]


def get_automorphic_structures(
        graph,
        atoms: tuple[int, ...],
    ) -> set[tuple[int, ...]]:
        """
        Returns all automorphic reoccurrences of the given atoms in structure
        Can be used for single atoms, two atoms (for bonds),
        three atoms(for angles), four atoms (for dihedrals) and n atoms for any
        structure for example identical reactive sides.

        :param atoms: Atoms to be used for search
        :type atoms: tuple[int, ...]
        :return: Automorphically equivalent sets of atoms
        :rtype: set[tuple[atom_id]]
        """
        return set(
            (
                tuple(mapping[atom] for atom in atoms)
                for mapping in graph.get_isomorphic_mappings(graph)
            )
        )