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
from stereomolgraph.algorithms.circular import color_refine_mg, color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    AtropBond,
    BondStereo,
    PlanarBond,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping


SiteKey = tuple[Literal["atom", "bond"], int | Bond]
StereoDescriptor = AtomStereo | BondStereo
SiteClass = Literal["specified", "distinct_colors", "unresolved", "geometry_only"]


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


def _candidate_descriptors(graph: StereoMolGraph) -> tuple[StereoDescriptor, ...]:
    """Use ONLY the stereo descriptors explicitly supplied by the input graph."""
    return (*graph.atom_stereo.values(), *graph.bond_stereo.values())


def detect_stereo_elements(graph: StereoMolGraph) -> StereoAnalysis:
    """Collect possible stereo sites and estimate stereo dependencies.

    Do not discard sites with equal ligand colors: another stereochemical
    assignment could make them stereogenic (para-stereogenicity).

    The input StereoMolGraph already contains ALL candidate geometry units.
    No stereo descriptor is inferred, added, or deleted during analysis.
    """
    topology = MolGraph(graph)
    base = _colors(topology)
    sites: list[StereoSite] = []

    for descriptor in _candidate_descriptors(graph):
        if descriptor.parity is None:
            states = tuple(
                sorted(
                    descriptor.get_isomers(),
                    key=lambda s: (type(s).__name__, repr(s.atoms), repr(s.parity)),
                )
            )
            classification: SiteClass = (
                "geometry_only"
                if len(states) == 1 and states[0].parity is None
                else (
                    "distinct_colors"
                    if _distinct_ligand_colors(_ligand_groups(descriptor), base)
                    else "unresolved"
                )
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
    unique: bool = True,
    max_isomers: int | None = None,
) -> Iterator[StereoMolGraph]:
    """Lazily enumerate local-state assignments; deduplicate at leaves only.

    No automorphism pruning is used. Outputs are candidate stereo graphs,
    not certified chemically feasible stereoisomers: unresolved sites
    could be inactive, and ring-geometry feasibility is not checked.
    Only existing descriptor slots are overwritten; no new sites are added.
    unique=True uses the existing StereoMolGraph graph isomorphism test.
    """
    if max_isomers is not None and max_isomers < 1:
        raise ValueError("max_isomers must be positive or None")
    if analysis is None:
        analysis = detect_stereo_elements(graph)

    # Reusing analysis from an unrelated graph must never insert new stereo units.
    input_descriptors = {
        _site_key(s): s for s in _candidate_descriptors(graph)
    }
    analysis_keys = [site.key for site in analysis.sites]
    if (
        len(analysis_keys) != len(input_descriptors)
        or set(analysis_keys) != set(input_descriptors)
    ):
        raise ValueError("analysis sites must match the input stereo descriptors")
    for site in analysis.sites:
        original = input_descriptors[site.key]
        if (
            type(site.template) is not type(original)
            or site.template.atoms != original.atoms
            or site.template.parity != original.parity
        ):
            raise ValueError("analysis descriptors must match the input graph")

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


def _replace_stereo_state(graph: StereoMolGraph, state: StereoDescriptor) -> None:
    """Overwrite the existing site descriptor with a local configuration."""
    key = _site_key(state)
    if key[0] == "atom":
        if key[1] not in graph.atom_stereo:
            raise ValueError("Cannot add a new atom stereo descriptor")
        graph.set_atom_stereo(state)
    else:
        if key[1] not in graph.bond_stereo:
            raise ValueError("Cannot add a new bond stereo descriptor")
        graph.set_bond_stereo(state)


def _anchored_stereo_isomorphic(
    first: StereoMolGraph, second: StereoMolGraph, key: SiteKey
) -> bool:
    """Exact stereo isomorphism requiring the designated site to map to itself.

    An atom site is individualized at its center; a bond site is individualized
    at both endpoints (which may still exchange with each other). This is an
    isomorphism test on complete stereo graphs, NOT McKay-style pruning.
    """
    if key[0] == "atom":
        anchors = frozenset((key[1],))
    else:
        anchors = key[1]

    def labels(g: StereoMolGraph) -> np.ndarray:
        atoms = tuple(g.atoms)
        seeds = [(int(g.get_atom_type(a)), a in anchors) for a in atoms]
        ordered = {seed: i for i, seed in enumerate(sorted(set(seeds)))}
        array = np.array([ordered[x] for x in seeds], dtype=np.int64)
        return color_refine_smg(g, atom_labels=array)

    return any(
        vf2pp_all_isomorphisms(
            first,
            second,
            stereo=True,
            atom_labels=(labels(first), labels(second)),
        )
    )


def locally_stereogenic_at_assignment(
    graph: StereoMolGraph,
    key: SiteKey,
) -> bool:
    """Determine whether a site has distinguishable local states at a fixed assignment.

    All *other* sites with more than one local state must be specified first.
    A site's alternative states are compared using an isomorphism which
    fixes the site, while keeping all other descriptors unchanged.

    This is exact for the *represented stereo graphs* assuming the underlying
    isomorphism implementation is correct. It is NOT a physical stability
    or geometry feasibility test, and it does not certify global conditional
    stereogenicity before other sites are assigned.
    """
    descriptors = _candidate_descriptors(graph)
    available = {_site_key(s): s for s in descriptors}
    if key not in available:
        raise ValueError("Site is not present among input stereo descriptors")

    for other in descriptors:
        if _site_key(other) == key:
            continue
        if other.parity is None and len(other.get_isomers()) > 1:
            raise ValueError(
                "All other variable stereo sites must be assigned first"
            )

    states = tuple(available[key].get_isomers())
    if len(states) < 2:
        return False

    # Fast sufficient result: strictly distinct constitutional ligand colors
    # imply that alternative local arrangements cannot be isomorphic at this site.
    base_colors = _colors(MolGraph(graph))
    if _distinct_ligand_colors(_ligand_groups(available[key]), base_colors):
        return True

    first = graph.copy()
    _replace_stereo_state(first, states[0])
    for state in states[1:]:
        alternate = graph.copy()
        _replace_stereo_state(alternate, state)
        if not _anchored_stereo_isomorphic(first, alternate, key):
            return True
    return False
