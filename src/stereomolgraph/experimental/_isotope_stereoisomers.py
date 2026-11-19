"""Automorphism-orbit-based isotope exchange and stereoisomer generation.

Atoms, including hydrogens, must be explicit. Exchange rules can specify
elements, alternative isotope mass numbers, and allowed atom IDs.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from numbers import Integral

import numpy as np

from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.algorithms.symmetry import atom_automorphism_classes
from stereomolgraph.graphs.mg import AtomId, MolGraph
from stereomolgraph.graphs.smg import StereoMolGraph
from stereomolgraph.periodic_table import PERIODIC_TABLE
from stereomolgraph.stereodescriptors import Tetrahedral


@dataclass(frozen=True)
class IsotopeExchange:
    """One exchange rule; positions refer to atom IDs, not array indices.

    atom_type: element symbol or atomic number, e.g. 'H' or 1.
    isotopes: possible target mass numbers, e.g. (2, 3) for D and T.
    source_isotope: original mass number, 0 means unspecified.
    positions: optional set of eligible atom IDs; None allows all matching.
    """

    atom_type: int | str
    isotopes: tuple[int, ...]
    source_isotope: int = 0
    positions: frozenset[AtomId] | None = None

    def __post_init__(self) -> None:
        try:
            atomic_number = int(PERIODIC_TABLE[self.atom_type])
        except (KeyError, TypeError) as exc:
            raise ValueError(f"Unknown atom type: {self.atom_type!r}") from exc
        object.__setattr__(self, "atom_type", atomic_number)
        source = MolGraph._validate_isotope(self.source_isotope)
        object.__setattr__(self, "source_isotope", source)
        try:
            values = tuple(self.isotopes)
        except TypeError as exc:
            raise ValueError("isotopes must be a nonempty sequence") from exc
        if not values:
            raise ValueError("isotopes must contain target mass numbers")
        targets = tuple(sorted(set(MolGraph._validate_isotope(x) for x in values)))
        if any(x == 0 or x == source for x in targets):
            raise ValueError("target isotopes must be positive and differ from source")
        object.__setattr__(self, "isotopes", targets)
        if self.positions is not None:
            positions = frozenset(self.positions)
            if any(not isinstance(a, Integral) or isinstance(a, bool) for a in positions):
                raise ValueError("positions must contain integer atom IDs")
            object.__setattr__(self, "positions", positions)


def _unique_graphs(graphs: Iterable[StereoMolGraph]) -> list[StereoMolGraph]:
    """Deduplicate completed stereo graphs by exact graph equality."""
    seen: set[StereoMolGraph] = set()
    unique = []
    for graph in graphs:
        graph.freeze()
        if graph not in seen:
            seen.add(graph)
            unique.append(graph.copy())
    return unique


def _exchange_options(
    graph: StereoMolGraph, rules: tuple[IsotopeExchange, ...]
) -> dict[AtomId, tuple[int, ...]]:
    """Precompute original site eligibility and all permitted target masses."""
    options: dict[AtomId, set[int]] = {atom: set() for atom in graph.atoms}
    atoms = set(graph.atoms)
    for rule in rules:
        if rule.positions is not None and not rule.positions <= atoms:
            raise ValueError(
                f"Unknown exchange atom IDs: {sorted(rule.positions - atoms)}"
            )
        for atom in graph.atoms:
            if rule.positions is not None and atom not in rule.positions:
                continue
            if int(graph.get_atom_type(atom)) != rule.atom_type:
                continue
            if graph.get_atom_isotope(atom) == rule.source_isotope:
                options[atom].update(rule.isotopes)
    return {atom: tuple(sorted(targets)) for atom, targets in options.items()}


class _OrbitLabels:
    """Stable integer labels for (element, current isotope, site options).

    Including the static per-atom allowed target masses prevents an
    automorphism from moving a partially labeled pattern to a position
    with different future substitution possibilities.
    """

    def __init__(
        self, graph: StereoMolGraph, options: dict[AtomId, tuple[int, ...]]
    ) -> None:
        self.options = options
        identities = {
            (int(graph.get_atom_type(atom)), mass, options[atom])
            for atom in graph.atoms
            for mass in (graph.get_atom_isotope(atom), *options[atom])
        }
        self.lookup = {item: index for index, item in enumerate(sorted(identities))}

    def signatures(self, graph: StereoMolGraph) -> tuple:
        return tuple(
            (int(graph.get_atom_type(atom)), graph.get_atom_isotope(atom),
             self.options[atom])
            for atom in graph.atoms
        )

    def labels(self, graph: StereoMolGraph) -> np.ndarray:
        values = self.signatures(graph)
        return np.fromiter(
            (self.lookup[item] for item in values),
            dtype=np.int64, count=len(values)
        )

    def invariant(self, graph: StereoMolGraph) -> tuple:
        return tuple(sorted(Counter(self.signatures(graph)).items()))


def _unique_patterns(
    candidates: Iterable[StereoMolGraph], labels: _OrbitLabels
) -> list[StereoMolGraph]:
    """Deduplicate constrained patterns across different parent paths."""
    buckets: dict[tuple, list[tuple[StereoMolGraph, np.ndarray]]] = defaultdict(list)
    unique = []
    for candidate in candidates:
        key = labels.invariant(candidate)
        candidate_colors = labels.labels(candidate)
        if any(
            next(
                vf2pp_all_isomorphisms(
                    candidate, previous, stereo=True,
                    atom_labels=(candidate_colors, previous_colors)
                ),
                None
            ) is not None
            for previous, previous_colors in buckets[key]
        ):
            continue
        buckets[key].append((candidate, candidate_colors))
        unique.append(candidate)
    return unique


def _orbit_representatives(
    graph: StereoMolGraph,
    originals: dict[AtomId, int],
    labels: _OrbitLabels,
) -> Iterator[AtomId]:
    """Return one eligible site from each automorphism orbit at this step."""
    colors = labels.labels(graph)
    mappings = vf2pp_all_isomorphisms(
        graph, graph, stereo=True, atom_labels=(colors, colors)
    )
    orbits = atom_automorphism_classes(graph, mappings=mappings)
    visited: set[AtomId] = set()
    for atom in sorted(graph.atoms):
        if (
            atom in visited
            or not labels.options[atom]
            or graph.get_atom_isotope(atom) != originals[atom]
        ):
            continue
        visited.update(orbits[atom])
        yield atom


def _isotope_patterns(
    graph: StereoMolGraph,
    options: dict[AtomId, tuple[int, ...]],
    min_substitutions: int,
    max_substitutions: int,
) -> Iterator[StereoMolGraph]:
    """At every depth compute new atom orbits and expand unique patterns."""
    original = {atom: graph.get_atom_isotope(atom) for atom in graph.atoms}
    labels = _OrbitLabels(graph, options)
    frontier = [graph.copy()]
    if min_substitutions == 0:
        yield graph.copy()

    for depth in range(1, max_substitutions + 1):
        def children() -> Iterator[StereoMolGraph]:
            for parent in frontier:
                for atom in _orbit_representatives(parent, original, labels):
                    for isotope in options[atom]:
                        child = parent.copy()
                        child.set_atom_attribute(atom, "isotope", isotope)
                        yield child

        frontier = _unique_patterns(children(), labels)
        if not frontier:
            return
        if depth >= min_substitutions:
            yield from frontier


def _equivalent_terminal_neighbors(graph: StereoMolGraph, atom: AtomId) -> bool:
    """Whether an exact terminal-atom swap makes tetrahedral stereo redundant.

    This is a conservative optimization. It does not assume that equal
    circular colors imply equivalent whole-molecule substituents.
    Other stereo descriptors referencing the terminal sites prevent pruning.
    """
    seen: set[tuple[int, int]] = set()
    for neighbor in graph.bonded_to(atom):
        if len(graph.bonded_to(neighbor)) != 1:
            continue
        if any(
            neighbor in descriptor.atoms
            for descriptor in graph.stereo.values()
            if not (isinstance(descriptor, Tetrahedral)
                    and descriptor.central_atom == atom)
        ):
            continue
        signature = (
            int(graph.get_atom_type(neighbor)),
            graph.get_atom_isotope(neighbor),
        )
        if signature in seen:
            return True
        seen.add(signature)
    return False


def _stereo_variables(
    graph: StereoMolGraph,
    infer_tetrahedral: bool,
) -> list[tuple[str, AtomId | frozenset[AtomId]]]:
    """Find undefined stereochemical elements; add conservative tetrahedral sites."""
    variables: list[tuple[str, AtomId | frozenset[AtomId]]] = []
    for atom in sorted(graph.atoms):
        descriptor = graph.get_atom_stereo(atom)
        if descriptor is None:
            if (
                infer_tetrahedral
                and int(graph.get_atom_type(atom)) in (6, 7, 14, 15, 32)
                and len(graph.bonded_to(atom)) == 4
            ):
                if _equivalent_terminal_neighbors(graph, atom):
                    continue
                graph.set_atom_stereo(
                    Tetrahedral((atom, *sorted(graph.bonded_to(atom))), None)
                )
                variables.append(("atom", atom))
        elif descriptor.parity is None:
            if isinstance(descriptor, Tetrahedral) and _equivalent_terminal_neighbors(
                graph, atom
            ):
                graph.delete_atom_stereo(atom)
            else:
                variables.append(("atom", atom))

    for bond, descriptor in sorted(
        graph.bond_stereo.items(), key=lambda pair: tuple(sorted(pair[0]))
    ):
        if descriptor.parity is None and len(descriptor.get_isomers()) > 1:
            variables.append(("bond", bond))
    return variables


def _get_descriptor(
    graph: StereoMolGraph, kind: str, key: AtomId | frozenset[AtomId]
):
    if kind == "atom":
        assert isinstance(key, int)
        return graph.get_atom_stereo(key)
    assert isinstance(key, frozenset)
    return graph.get_bond_stereo(key)


def _set_descriptor(graph: StereoMolGraph, kind: str, descriptor) -> None:
    if kind == "atom":
        graph.set_atom_stereo(descriptor)
    else:
        graph.set_bond_stereo(descriptor)


def _remove_redundant_stereo(
    graph: StereoMolGraph,
    variables: list[tuple[str, AtomId | frozenset[AtomId]]],
) -> StereoMolGraph:
    """Remove descriptors whose two configurations are isomorphic.

    Only introduced or originally undefined descriptors may be removed.
    Pre-existing, assigned stereochemistry is always conserved.
    """
    result = graph.copy()
    changed = True
    while changed:
        changed = False
        for kind, key in variables:
            descriptor = _get_descriptor(result, kind, key)
            if descriptor is None:
                continue
            alternatives = [
                option
                for option in descriptor.get_isomers()
                if option != descriptor
            ]
            if len(alternatives) != 1:
                continue
            alternative = result.copy()
            _set_descriptor(alternative, kind, alternatives[0])
            if alternative == result:
                if kind == "atom":
                    assert isinstance(key, int)
                    result.delete_atom_stereo(key)
                else:
                    assert isinstance(key, frozenset)
                    result.delete_bond_stereo(key)
                changed = True
    return result


def _stereo_variants(
    isotope_graph: StereoMolGraph, infer_tetrahedral: bool
) -> Iterator[StereoMolGraph]:
    """Generate full assignments, then deduplicate exact stereoisomers.

    Do not deduplicate intermediate stereo assignments: undefined parity can
    temporarily hide dependent stereogenicity, and early pruning can lose
    valid stereoisomers.
    """
    template = isotope_graph.copy()
    variables = _stereo_variables(template, infer_tetrahedral)
    options: list[tuple] = []
    for kind, key in variables:
        descriptor = _get_descriptor(template, kind, key)
        assert descriptor is not None
        options.append(tuple(sorted(descriptor.get_isomers(), key=repr)))

    def candidates() -> Iterator[StereoMolGraph]:
        for assignment in itertools.product(*options):
            result = template.copy()
            for (kind, _key), descriptor in zip(variables, assignment):
                _set_descriptor(result, kind, descriptor)
            yield _remove_redundant_stereo(result, variables)

    yield from _unique_graphs(candidates())


def generate_isotope_stereoisomers(
    graph: StereoMolGraph,
    exchanges: IsotopeExchange | Iterable[IsotopeExchange],
    *,
    min_substitutions: int = 1,
    max_substitutions: int | None = None,
    infer_tetrahedral: bool = True,
) -> Iterator[StereoMolGraph]:
    """Generate each unique isotopically substituted stereoisomer exactly once.

    One or more IsotopeExchange rules give the element, target isotope choices,
    source isotope, and optionally the allowed atom IDs. Each atom is
    exchanged at most once. The substitution limits count modified atoms
    over all rules. Use explicit atoms for isotope exchange (including H).

    At each substitution depth, automorphisms of the partially substituted
    graph give atom orbits. Orbit labels include allowed target options to
    preserve positional restrictions. A second exact isomorphism filter
    removes duplicates reached through different augmentation orders.

    Undefined stereochemical descriptors and newly inferred tetrahedral
    elements are enumerated without inverting originally defined parity.
    A final global exact stereo-aware set removes duplicates even between
    different isotope rules. Both distinct enantiomers are returned.
    """
    if not isinstance(graph, StereoMolGraph):
        raise TypeError("graph must be a StereoMolGraph")
    rules = (exchanges,) if isinstance(exchanges, IsotopeExchange) else tuple(exchanges)
    if not rules or any(not isinstance(rule, IsotopeExchange) for rule in rules):
        raise ValueError("exchanges must contain one or more IsotopeExchange rules")
    if (
        not isinstance(min_substitutions, int)
        or isinstance(min_substitutions, bool)
        or min_substitutions < 0
    ):
        raise ValueError("min_substitutions must be a non-negative integer")

    options = _exchange_options(graph, rules)
    eligible_count = sum(bool(targets) for targets in options.values())
    if max_substitutions is None:
        max_substitutions = eligible_count
    if (
        not isinstance(max_substitutions, int)
        or isinstance(max_substitutions, bool)
        or not min_substitutions <= max_substitutions <= eligible_count
    ):
        raise ValueError(
            "Expected min_substitutions <= max_substitutions <= eligible atom count"
        )

    seen: set[StereoMolGraph] = set()
    for pattern in _isotope_patterns(
        graph, options, min_substitutions, max_substitutions
    ):
        for product in _stereo_variants(pattern, infer_tetrahedral):
            product.freeze()
            if product not in seen:
                seen.add(product)
                yield product.copy()
