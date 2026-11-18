"""Symmetry-distinct isotope exchange and stereoisomer generation.

Explicit atoms are required for isotope exchange (including hydrogens).
The stereo enumerator uses all undefined descriptors supplied by the caller
and can infer tetrahedral sites for ordinary four-coordinate main-group atoms.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Iterator

from stereomolgraph.graphs.mg import AtomId
from stereomolgraph.graphs.smg import StereoMolGraph
from stereomolgraph.periodic_table import PERIODIC_TABLE
from stereomolgraph.stereodescriptors import Tetrahedral


def _unique_graphs(graphs: Iterable[StereoMolGraph]) -> list[StereoMolGraph]:
    """Deduplicate using the library's hash and exact stereo isomorphism."""
    seen: set[StereoMolGraph] = set()
    unique = []
    for graph in graphs:
        graph.freeze()
        if graph not in seen:
            seen.add(graph)
            unique.append(graph.copy())
    return unique


def _isotope_patterns(
    graph: StereoMolGraph,
    sites: tuple[AtomId, ...],
    source_isotope: int,
    target_isotope: int,
    min_substitutions: int,
    max_substitutions: int,
) -> Iterator[StereoMolGraph]:
    """Generate isotope-placement orbit representatives at each depth.

    Each frontier contains only unique isotope patterns. Exact graph
    isomorphism collapses all substitution orders and symmetry-equivalent
    sites, without enumerating the potentially factorial automorphism group.
    """
    frontier = [graph.copy()]
    if min_substitutions == 0:
        yield graph.copy()

    for count in range(1, max_substitutions + 1):
        def augment() -> Iterator[StereoMolGraph]:
            for parent in frontier:
                for atom in sites:
                    if parent.get_atom_isotope(atom) != source_isotope:
                        continue
                    child = parent.copy()
                    child.set_atom_attribute(atom, "isotope", target_isotope)
                    yield child

        frontier = _unique_graphs(augment())
        if not frontier:
            return
        if count >= min_substitutions:
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
    element: int | str,
    isotope: int,
    *,
    source_isotope: int = 0,
    min_substitutions: int = 1,
    max_substitutions: int | None = None,
    atoms: Iterable[AtomId] | None = None,
    infer_tetrahedral: bool = True,
) -> Iterator[StereoMolGraph]:
    """Enumerate unique isotope-exchanged stereoisomers.

    Example: H -> D at one or more explicit H atoms::

        products = generate_isotope_stereoisomers(
            graph, "H", 2, min_substitutions=1, max_substitutions=2
        )

    Isotope=0 denotes an unspecified source isotope; to replace an explicitly
    labeled isotope, pass its mass number as `source_isotope`.

    Only atoms matching the requested element and source isotope are exchanged.
    Stereochemical constraints with assigned parity are preserved. Undefined
    stereo descriptors are enumerated, and missing tetrahedral descriptors are
    inferred for four-coordinate C, N, Si, P, and Ge atoms when enabled.
    Other stereo element types must have appropriate descriptors in the input.

    Returns each distinct molecule once, including both enantiomers. Requires
    explicit exchangeable atoms; it does not add implicit hydrogens or infer
    other geometries from bond orders alone. Isotope patterns and resulting
    stereo configurations are compared by exact graph isomorphism after hash
    filtering. The number of distinct products can grow exponentially.
    """
    if not isinstance(graph, StereoMolGraph):
        raise TypeError("graph must be a StereoMolGraph")
    atomic_number = int(PERIODIC_TABLE[element])
    source_isotope = graph._validate_isotope(source_isotope)
    isotope = graph._validate_isotope(isotope)
    if not isotope or isotope == source_isotope:
        raise ValueError("target isotope must be positive and differ from source")
    if not isinstance(min_substitutions, int) or isinstance(min_substitutions, bool):
        raise ValueError("min_substitutions must be an integer")
    if min_substitutions < 0:
        raise ValueError("min_substitutions must be non-negative")

    allowed = None if atoms is None else set(atoms)
    if allowed is not None and not allowed <= set(graph.atoms):
        raise ValueError("atoms must be a subset of graph.atoms")

    sites = tuple(
        sorted(
            atom for atom in graph.atoms
            if (allowed is None or atom in allowed)
            and int(graph.get_atom_type(atom)) == atomic_number
            and graph.get_atom_isotope(atom) == source_isotope
        )
    )
    if max_substitutions is None:
        max_substitutions = len(sites)
    if (
        not isinstance(max_substitutions, int)
        or isinstance(max_substitutions, bool)
        or max_substitutions < min_substitutions
        or max_substitutions > len(sites)
    ):
        raise ValueError(
            "max_substitutions must be an integer between "
            "min_substitutions and the number of eligible atoms"
        )

    for pattern in _isotope_patterns(
        graph, sites, source_isotope, isotope,
        min_substitutions, max_substitutions
    ):
        yield from _stereo_variants(pattern, infer_tetrahedral)
