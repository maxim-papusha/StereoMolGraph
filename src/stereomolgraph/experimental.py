from __future__ import annotations

import itertools
from collections import deque
from typing import TYPE_CHECKING

from stereomolgraph import StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.graphs.scrg import Change

if TYPE_CHECKING:
    from collections.abc import Collection, Iterable
    from typing import Optional

    from stereomolgraph.graphs.smg import AtomId, Bond, StereoMolGraph


def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: Optional[Iterable[AtomId]] = None,
    bonds: Optional[Iterable[Bond]] = None,
) -> Collection[StereoMolGraph]:
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

    isomers: set[StereoMolGraph] = set()
    enantiomers_set: set[StereoMolGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)

        if stereoisomer not in enantiomers_set:
            isomers.add(stereoisomer)

            if not enantiomers:
                enantiomers_set.add(stereoisomer.enantiomer())

    return isomers


def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: Optional[Iterable[AtomId]] = None,
    bonds: Optional[Iterable[Bond]] = None,
) -> Collection[StereoCondensedReactionGraph]:
    # TODO: extend to more than fleeting stereochemistry
    # add checks if fleeting stereochemistry is valid relative to formed
    # and broken

    if atoms is None:
        atom_stereos = [
            stereo.get_isomers()
            for a in graph.atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        ]
    else:
        atom_stereos = [
            stereo.get_isomers()
            for a in atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        ]

    if bonds is None:
        bond_stereos = [
            stereo.get_isomers()
            for b in graph.bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        ]
    else:
        bond_stereos = [
            stereo.get_isomers()
            for b in bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        ]

    isomers: list[StereoCondensedReactionGraph] = []
    enantiomers_set: set[StereoCondensedReactionGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)

        if stereoisomer not in enantiomers_set:
            isomers.append(stereoisomer)

            if not enantiomers:
                enantiomers_set.add(stereoisomer.enantiomer())

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
