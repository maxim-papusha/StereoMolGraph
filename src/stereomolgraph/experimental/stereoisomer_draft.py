"""Experimental stereoisomer candidate detection and enumeration.

This is an exploratory algorithm, not a complete stereogenicity or chemical
feasibility validator. It does NOT prune using graph automorphisms.
See docs/stereoisomer-algorithm-draft.md for scope and limitations.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from stereomolgraph import Bond, MolGraph, StereoMolGraph
from stereomolgraph.algorithms.circular import color_refine_mg
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    AtropBond,
    BondStereo,
    PlanarBond,
    Tetrahedral,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping


SiteKey = tuple[Literal["atom", "bond"], int | Bond]
StereoDescriptor = AtomStereo | BondStereo
SiteClass = Literal["specified", "distinct_colors", "unresolved"]


@dataclass(frozen=True)
class StereoSite:
    """One geometry-qualified stereo candidate with its local state domain."""

    key: SiteKey
    template: StereoDescriptor
    states: tuple[StereoDescriptor, ...]
    classification: SiteClass


@dataclass(frozen=True)
class StereoAnalysis:
    """Candidate sites and heuristic directed stereo-influence edges.

    The dependencies[target] set contains sites whose positional marking
    splits an initially equal-color ligand pair at target. These are only
    diagnostics, not a certified dependency graph.
    """

    sites: tuple[StereoSite, ...]
    dependencies: Mapping[SiteKey, frozenset[SiteKey]]


def _site_key(stereo: StereoDescriptor) -> SiteKey:
    if isinstance(stereo, AtomStereo):
        return ("atom", stereo.central_atom)
    if isinstance(stereo, BondStereo):
        return ("bond", Bond(stereo.bond))
    raise TypeError(f"Unsupported descriptor: {type(stereo)!r}")


def _ligand_groups(stereo: StereoDescriptor) -> tuple[tuple[int | None, ...], ...]:
    """Ligands compared independently at each end of a stereo unit."""
    if isinstance(stereo, (PlanarBond, AtropBond)):
        return (tuple(stereo.atoms[:2]), tuple(stereo.atoms[4:6]))
    if isinstance(stereo, AtomStereo):
        return (tuple(stereo.atoms[1:]),)
    return ()  # Rigid-bond classes are not subject to this quick check.


def _colors(graph: MolGraph, marked: frozenset[int] = frozenset()) -> dict[int, int]:
    """1-WL refinement from atomic numbers, individualizing a marked site."""
    atoms = tuple(graph.atoms)
    seed = [(int(graph.get_atom_type(a)), int(a in marked)) for a in atoms]
    rank = {value: i for i, value in enumerate(sorted(set(seed)))}
    initial = np.array([rank[x] for x in seed], dtype=np.int64)
    refined = color_refine_mg(graph, atom_labels=initial)
    return {a: int(c) for a, c in zip(atoms, refined, strict=True)}


def _pair_partition_breaks(
    groups: tuple[tuple[int | None, ...], ...],
    initial: Mapping[int, int],
    marked: Mapping[int, int],
) -> bool:
    def color(atom: int | None, palette: Mapping[int, int]) -> tuple[str, int]:
        return ("missing", 0) if atom is None else ("atom", palette[atom])

    for group in groups:
        for i, j in itertools.combinations(group, 2):
            if color(i, initial) == color(j, initial) and (
                color(i, marked) != color(j, marked)
            ):
                return True
    return False


def _distinct_ligand_colors(
    groups: tuple[tuple[int | None, ...], ...],
    palette: Mapping[int, int],
) -> bool:
    if not groups:
        return False
    for group in groups:
        signatures = [
            ("missing", 0) if a is None else ("atom", palette[a])
            for a in group
        ]
        if len(signatures) != len(set(signatures)):
            return False
    return True


def _candidate_descriptors(
    graph: StereoMolGraph, *, infer_tetrahedral: bool
) -> tuple[StereoDescriptor, ...]:
    candidates: list[StereoDescriptor] = [
        *graph.atom_stereo.values(), *graph.bond_stereo.values()
    ]
    if infer_tetrahedral:
        for a in graph.atoms:
            # Conservative structural shortcut: only explicitly bonded,
            # four-coordinate C/Si. Other geometries need a template.
            if (
                a not in graph.atom_stereo
                and graph.get_atom_type(a) in (6, 14)
                and len(graph.bonded_to(a)) == 4
            ):
                candidates.append(
                    Tetrahedral((a, *sorted(graph.bonded_to(a))), None)
                )
    return tuple(candidates)


def detect_stereo_elements(
    graph: StereoMolGraph,
    *,
    infer_tetrahedral: bool = False,
) -> StereoAnalysis:
    """Collect possible stereo sites and estimate stereo dependencies.

    Do not discard sites with equal ligand colors: another stereochemical
    assignment could make them stereogenic (para-stereogenicity).

    Normally geometry candidates are already present on the graph.
    infer_tetrahedral=True also nominates four-coordinate C/Si sites; it is
    not general chemical stereoperception.
    """
    topology = MolGraph(graph)
    base = _colors(topology)
    sites: list[StereoSite] = []

    for descriptor in _candidate_descriptors(
        graph, infer_tetrahedral=infer_tetrahedral
    ):
        if descriptor.parity is None:
            states = tuple(
                sorted(
                    descriptor.get_isomers(),
                    key=lambda s: repr((type(s).__name__, s.canonical_form())),
                )
            )
            classification: SiteClass = (
                "distinct_colors"
                if _distinct_ligand_colors(_ligand_groups(descriptor), base)
                else "unresolved"
            )
        else:
            states = (descriptor,)
            classification = "specified"
        if not states:
            raise ValueError(f"No local states for {descriptor!r}")
        sites.append(StereoSite(_site_key(descriptor), descriptor, states, classification))

    influence: dict[SiteKey, set[SiteKey]] = {site.key: set() for site in sites}
    for source in sites:
        site_atoms = (
            frozenset((source.key[1],))
            if source.key[0] == "atom"
            else source.key[1]
        )
        marked = _colors(topology, site_atoms)
        for target in sites:
            if source.key == target.key or target.classification != "unresolved":
                continue
            if _pair_partition_breaks(
                _ligand_groups(target.template), base, marked
            ):
                influence[target.key].add(source.key)

    return StereoAnalysis(
        tuple(sites),
        {key: frozenset(values) for key, values in influence.items()},
    )


def enumerate_stereo_assignments(
    graph: StereoMolGraph,
    analysis: StereoAnalysis | None = None,
    *,
    infer_tetrahedral: bool = False,
    unique: bool = True,
    max_isomers: int | None = None,
) -> Iterator[StereoMolGraph]:
    """Lazily enumerate local-state assignments; deduplicate at leaves only.

    No automorphism pruning is used. Outputs are candidate stereo graphs,
    not certified chemically feasible stereoisomers: unresolved sites
    could be inactive, and ring-geometry feasibility is not checked.
    unique=True uses the existing StereoMolGraph graph isomorphism test.
    """
    if max_isomers is not None and max_isomers < 1:
        raise ValueError("max_isomers must be positive or None")
    if analysis is None:
        analysis = detect_stereo_elements(
            graph, infer_tetrahedral=infer_tetrahedral
        )

    seen: set[StereoMolGraph] = set()
    count = 0
    for assignments in itertools.product(*(s.states for s in analysis.sites)):
        candidate = graph.copy()
        for stereo in assignments:
            if isinstance(stereo, AtomStereo):
                candidate.set_atom_stereo(stereo)
            elif isinstance(stereo, BondStereo):
                candidate.set_bond_stereo(stereo)
            else:
                raise TypeError(f"Unsupported descriptor: {type(stereo)!r}")

        if unique:
            candidate.freeze()
            if candidate in seen:
                continue
            seen.add(candidate)
            candidate = candidate.copy(frozen=False)

        yield candidate
        count += 1
        if max_isomers is not None and count >= max_isomers:
            return
