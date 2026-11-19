"""Isotope exchange with symmetry-aware stereoisomer enumeration.

Exchangeable atoms, including hydrogen, must be explicit. Only fully
constructed isotope patterns are returned; all duplicates are removed.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from itertools import product
from numbers import Integral

import numpy as np

from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.algorithms.symmetry import atom_automorphism_classes
from stereomolgraph.graphs.mg import AtomId, MolGraph
from stereomolgraph.graphs.smg import StereoMolGraph
from stereomolgraph.periodic_table import PERIODIC_TABLE
from stereomolgraph.stereodescriptors import Tetrahedral


@dataclass(frozen=True, slots=True)
class IsotopeExchange:
    """Allowed isotope replacements for an element.

    :param atom_type: Element symbol or atomic number, e.g. "H" or 1.
    :param isotopes: Possible target mass numbers, e.g. (2, 3) for D and T.
    :param source_isotope: Original isotope (0 means unspecified, not 1H).
    :param positions: Optional allowed atom IDs; None means every matching atom.

    A site may be exchanged only once. Multiple rules may overlap: their
    target isotope options are combined at each eligible atom.
    """

    atom_type: int | str
    isotopes: tuple[int, ...]
    source_isotope: int = 0
    positions: frozenset[AtomId] | None = None

    def __post_init__(self) -> None:
        try:
            element = PERIODIC_TABLE[self.atom_type]
        except (KeyError, TypeError) as exc:
            raise ValueError(f"Unknown atom type: {self.atom_type!r}") from exc

        source = MolGraph._validate_isotope(self.source_isotope)
        try:
            targets = tuple(sorted({
                MolGraph._validate_isotope(mass) for mass in self.isotopes
            }))
        except TypeError as exc:
            raise ValueError("isotopes must be a nonempty iterable") from exc
        if not targets or any(mass == 0 or mass == source for mass in targets):
            raise ValueError("Target isotopes must be positive and differ from source")

        positions = None if self.positions is None else frozenset(self.positions)
        if positions is not None and any(
            not isinstance(atom, Integral) or isinstance(atom, bool)
            for atom in positions
        ):
            raise ValueError("positions must contain integer atom IDs")

        object.__setattr__(self, "atom_type", int(element))
        object.__setattr__(self, "source_isotope", source)
        object.__setattr__(self, "isotopes", targets)
        object.__setattr__(self, "positions", positions)


def _isotope_patterns(
    graph: StereoMolGraph,
    allowed: dict[AtomId, tuple[int, ...]],
    substitutions: int,
) -> Iterator[StereoMolGraph]:
    """Construct patterns using atom orbits at every substitution step.

    Atom labels also encode the permitted targets at each position. These
    restrictions must be preserved while deduplicating partial patterns:
    two chemically identical partial patterns may have different remaining
    substitution options.
    """
    initial = {atom: graph.get_atom_isotope(atom) for atom in graph.atoms}
    elements = {atom: int(graph.get_atom_type(atom)) for atom in graph.atoms}
    identities = {
        (elements[atom], mass, allowed.get(atom, ()))
        for atom in graph.atoms
        for mass in (initial[atom], *allowed.get(atom, ()))
    }
    color_of = {
        identity: index for index, identity in enumerate(sorted(identities))
    }

    def colors(molecule: StereoMolGraph) -> np.ndarray:
        return np.array([
            color_of[(
                elements[atom],
                molecule.get_atom_isotope(atom),
                allowed.get(atom, ()),
            )]
            for atom in molecule.atoms
        ], dtype=np.int64)

    frontier = [graph.copy()]
    for _ in range(substitutions):
        # Isomorphism buckets are reset for each substitution depth.
        buckets: dict[tuple[int, ...], list[tuple[StereoMolGraph, np.ndarray]]] = {}
        next_frontier: list[StereoMolGraph] = []

        for parent in frontier:
            parent_colors = colors(parent)
            automorphisms = vf2pp_all_isomorphisms(
                parent, parent, stereo=True,
                atom_labels=(parent_colors, parent_colors),
            )
            orbits = atom_automorphism_classes(parent, mappings=automorphisms)
            visited: set[AtomId] = set()

            for atom in sorted(allowed):
                if atom in visited or parent.get_atom_isotope(atom) != initial[atom]:
                    continue
                visited.update(orbits[atom])

                for mass in allowed[atom]:
                    child = parent.copy()
                    child.set_atom_attribute(atom, "isotope", mass)
                    child_colors = colors(child)
                    bucket = buckets.setdefault(
                        tuple(sorted(child_colors.tolist())), []
                    )
                    if any(
                        next(vf2pp_all_isomorphisms(
                            child, previous, stereo=True,
                            atom_labels=(child_colors, previous_colors),
                        ), None) is not None
                        for previous, previous_colors in bucket
                    ):
                        continue
                    bucket.append((child, child_colors))
                    next_frontier.append(child)

        frontier = next_frontier
        if not frontier:
            return

    yield from frontier


def _stereo_variants(
    graph: StereoMolGraph, infer_tetrahedral: bool
) -> Iterator[StereoMolGraph]:
    """Complete undefined stereochemistry and remove nonstereogenic descriptors.

    All stereo variables are assigned before testing equivalence; otherwise,
    dependent stereogenic elements can be inadvertently removed.
    """
    template = graph.copy()
    variables = []

    for atom in sorted(template.atoms):
        descriptor = template.get_atom_stereo(atom)
        if (
            descriptor is None
            and infer_tetrahedral
            and int(template.get_atom_type(atom)) in (6, 7, 14, 15, 32)
            and len(template.bonded_to(atom)) == 4
        ):
            neighbors = tuple(sorted(template.bonded_to(atom)))

            # Two identical terminal ligands cannot create tetrahedral stereo.
            terminal = [
                (int(template.get_atom_type(n)), template.get_atom_isotope(n))
                for n in neighbors
                if len(template.bonded_to(n)) == 1
                and not any(n in stereo.atoms for stereo in template.stereo.values())
            ]
            if len(terminal) == len(set(terminal)):
                descriptor = Tetrahedral((atom, *neighbors), None)
                template.set_atom_stereo(descriptor)

        if descriptor is not None and descriptor.parity is None:
            variables.append(("atom", atom, tuple(descriptor.get_isomers())))

    for bond, descriptor in sorted(
        template.bond_stereo.items(), key=lambda item: tuple(sorted(item[0]))
    ):
        if descriptor.parity is None:
            isomers = tuple(descriptor.get_isomers())
            if len(isomers) > 1:
                variables.append(("bond", bond, isomers))

    seen: set[StereoMolGraph] = set()
    for assignment in product(*(options for _, _, options in variables)):
        candidate = template.copy()
        for (kind, _key, _options), descriptor in zip(variables, assignment):
            if kind == "atom":
                candidate.set_atom_stereo(descriptor)
            else:
                candidate.set_bond_stereo(descriptor)

        # A descriptor with identical opposite configurations is redundant.
        # Repeat because removing one may make another descriptor redundant.
        changed = True
        while changed:
            changed = False
            for kind, key, _options in variables:
                descriptor = (
                    candidate.get_atom_stereo(key)
                    if kind == "atom"
                    else candidate.get_bond_stereo(key)
                )
                if descriptor is None:
                    continue
                opposites = [
                    option for option in descriptor.get_isomers()
                    if option != descriptor
                ]
                if len(opposites) != 1:
                    continue

                opposite = candidate.copy()
                if kind == "atom":
                    opposite.set_atom_stereo(opposites[0])
                else:
                    opposite.set_bond_stereo(opposites[0])

                if candidate == opposite:
                    if kind == "atom":
                        candidate.delete_atom_stereo(key)
                    else:
                        candidate.delete_bond_stereo(key)
                    changed = True

        candidate.freeze()
        if candidate not in seen:
            seen.add(candidate)
            yield candidate.copy()


def generate_isotope_stereoisomers(
    graph: StereoMolGraph,
    exchanges: IsotopeExchange | Iterable[IsotopeExchange],
    *,
    substitutions: int | None = None,
    infer_tetrahedral: bool = True,
) -> Iterator[StereoMolGraph]:
    """Yield unique completed isotope-exchange stereoisomers.

    Each rule specifies the element, permitted isotope(s) and optional atom
    IDs. By default, all eligible atoms are substituted. Use
    `substitutions=n` for exactly n exchanges, not 0 through n.

    Each construction step computes atom orbits of the partially labeled
    graph and extends one representative per orbit. Exact isomorphism removes
    equivalent substitution paths; a final stereo-aware uniqueness check
    ensures that no molecule is returned twice. Originally fixed stereo
    parities are preserved, and distinct enantiomers are retained.

    All exchangeable atoms (including H) must be explicit. Missing
    tetrahedral descriptors are inferred only for four-coordinate
    C/N/Si/P/Ge; other stereo types require input descriptors.
    """
    if not isinstance(graph, StereoMolGraph):
        raise TypeError("graph must be a StereoMolGraph")
    rules = (exchanges,) if isinstance(exchanges, IsotopeExchange) else tuple(exchanges)
    if not rules or any(not isinstance(rule, IsotopeExchange) for rule in rules):
        raise ValueError("exchanges must contain one or more IsotopeExchange rules")

    options: dict[AtomId, set[int]] = {atom: set() for atom in graph.atoms}
    for rule in rules:
        if rule.positions is not None and not rule.positions <= set(graph.atoms):
            raise ValueError(
                "Unknown exchange atom IDs: "
                f"{sorted(rule.positions - set(graph.atoms))}"
            )
        for atom in graph.atoms:
            if (
                (rule.positions is None or atom in rule.positions)
                and int(graph.get_atom_type(atom)) == rule.atom_type
                and graph.get_atom_isotope(atom) == rule.source_isotope
            ):
                options[atom].update(rule.isotopes)

    allowed = {
        atom: tuple(sorted(targets)) for atom, targets in options.items() if targets
    }
    if substitutions is None:
        # An empty exchange specification has no completed products.
        if not allowed:
            return
        substitutions = len(allowed)
    if (
        not isinstance(substitutions, int)
        or isinstance(substitutions, bool)
        or not 0 <= substitutions <= len(allowed)
    ):
        raise ValueError(
            f"substitutions must be an integer between 0 and {len(allowed)}"
        )

    seen: set[StereoMolGraph] = set()
    for pattern in _isotope_patterns(graph, allowed, substitutions):
        for product in _stereo_variants(pattern, infer_tetrahedral):
            product.freeze()
            if product not in seen:
                seen.add(product)
                yield product.copy()
