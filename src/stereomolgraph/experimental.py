from __future__ import annotations

import itertools
from collections import deque
from typing import TYPE_CHECKING

from stereomolgraph import (
    AtomId,
    Bond,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms.color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.graphs.scrg import Change


if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

def unique_generator(input_generator: Iterator) -> Iterator:
    """
    A generator that yields unique objects from another generator.
    
    Args:
        input_generator: A generator yielding hashable objects
        
    Yields:
        Only the first occurrence of each unique object from the input generator
    """
    seen_hash = set()
    for item in input_generator:
        item_hash = hash(item)
        if item_hash not in seen_hash:
            seen_hash.add(item_hash)
            yield item

def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoMolGraph]:
    """Generates all unique stereoisomers of a StereoMolGraph by generation of
    all combinations of parities. Only includes stereocenters which have a
    parity of None. If a parity is set, it is not changed.

    If include_enantiomers is True, both enantiomers of a stereoisomer are
    included, if it is False, only one enantiomer is included.

    Args:
        enantiomers: If True, both enantiomers are included
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism
        
    Yields:
        StereoMolGraph: Each unique stereoisomer (and enantiomer if requested)
    """
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if ((stereo := graph.get_atom_stereo(a)) and stereo.parity is None)
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (stereo := graph.get_atom_stereo(a)) is not None
        )

    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if ((stereo := graph.get_bond_stereo(b)) and stereo.parity is None)
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (stereo := graph.get_bond_stereo(b)) is not None
        )

    seen_hash: set[int] = set()
    enantiomers_seen_hash: set[int] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)

        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_hash:
            seen_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in enantiomers_seen_hash:
                    enantiomers_seen_hash.add(enantiomer_hash)
                    yield enantiomer

def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoCondensedReactionGraph]:
    """Generates all unique fleeting stereoisomers of a StereoCondensedReactionGraph.
    
    Only includes stereocenters which have a parity of None for the fleeting change.
    If a parity is set, it is not changed.

    Args:
        graph: The reaction graph to generate isomers from
        enantiomers: If True, both enantiomers are included (default: True)
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism
        
    Yields:
        StereoCondensedReactionGraph: Each unique fleeting stereoisomer
    """
    # Get atom stereoisomers
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    # Get bond stereoisomers
    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    seen_isomers_hash = set()
    seen_enantiomers_hash = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)
        
        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_isomers_hash:
            seen_isomers_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in seen_enantiomers_hash:
                    seen_enantiomers_hash.add(enantiomer_hash)
                    yield enantiomer


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
    colorings = color_refine_smg(graph)
    mappings = vf2pp_all_isomorphisms(graph, graph, atom_labels=(colorings, colorings), stereo=True)
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]
