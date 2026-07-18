from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from stereomolgraph import AtomId, Bond, StereoCondensedReactionGraph, StereoMolGraph
from stereomolgraph.graphs.crg import Change

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator


def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = False,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoMolGraph]:
    """
    Generate all unique stereoisomers of a StereoMolGraph.

    Stereoisomers are constructed by enumerating all combinations of
    stereochemical parities. Only stereocenters with undefined parity
    (``None``) are considered; predefined parities remain unchanged.

    If ``enantiomers`` is ``False``, both enantiomers of each stereoisomer
    are returned. Otherwise, only one representative per enantiomeric pair
    is included.

    :param enantiomers: Whether to remove enantiomers.
    :param atoms: Optional subset of atoms to consider for stereoisomerism.
    :param bonds: Optional subset of bonds to consider for stereoisomerism.

    :yields: Unique stereoisomers of the input graph.
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

    seen: set[StereoMolGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)
        stereoisomer.freeze()
        if stereoisomer in seen:
            continue

        enantiomer = stereoisomer.enantiomer()
        enantiomer.freeze()

        seen.add(stereoisomer)
        seen.add(enantiomer)

        yield stereoisomer.copy(frozen=False)

        # enantiomers=False means include both members of each enantiomeric pair.
        if not enantiomers and enantiomer != stereoisomer:
            yield enantiomer.copy(frozen=False)


def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoCondensedReactionGraph]:
    """
    Generate all unique fleeting stereoisomers of a
    StereoCondensedReactionGraph.

    Stereoisomers are constructed by enumerating all combinations of
    fleeting stereochemical parities. Only fleeting stereocenters with
    undefined parity (``None``) are considered; predefined parities remain
    unchanged.

    If ``enantiomers`` is ``False``, both enantiomers of each fleeting
    stereoisomer are returned. Otherwise, only one representative per
    enantiomeric pair is included.

    :param enantiomers: Whether to remove enantiomers.
    :param atoms: Optional subset of atoms to consider for stereoisomerism.
    :param bonds: Optional subset of bonds to consider for stereoisomerism.

    :yields: Unique fleeting stereoisomers of the input graph.
    """
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

    seen: set[StereoCondensedReactionGraph] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)

        stereoisomer.freeze()
        if stereoisomer in seen:
            continue

        enantiomer = stereoisomer.enantiomer()
        enantiomer.freeze()

        seen.add(stereoisomer)
        seen.add(enantiomer)

        yield stereoisomer

        if not enantiomers and enantiomer != stereoisomer:
            yield enantiomer
