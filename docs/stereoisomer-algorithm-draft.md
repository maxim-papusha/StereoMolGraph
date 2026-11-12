# Draft: stereogenicity detection and stereoisomer enumeration

**Status:** research prototype, not a chemically complete or formally validated
stereoisomer generator. This design is intentionally separate from the existing
enumerator. It contains **no automorphism-based enumeration pruning**.

Prototype: [`src/stereomolgraph/experimental/stereoisomer_draft.py`](../src/stereomolgraph/experimental/stereoisomer_draft.py).

## 1. Definitions and correctness criteria

There are four different notions. They should **not** be conflated:

1. **Geometrically eligible stereo unit:** a central atom, bond, or higher-order
   unit for which a local coordination geometry is chemically meaningful
   (tetrahedral, square planar, octahedral, planar double bond, atropisomer, etc.).
2. **Potentially stereogenic unit:** a unit that may distinguish configurations
   for at least one allowed assignment of all other units. This includes
   para-stereogenic and mutually dependent configurations.
3. **Currently stereogenic unit:** a unit with two local states distinguishable
   when the other stereochemical information is fixed.
4. **Distinct molecular stereoisomer:** an entire stereo-configured molecular
   graph inequivalent to any other under the chosen stereochemical isomorphism.

A fifth notion, **physical feasibility**, is not determined by the graph
alone (e.g. strained ring trans alkene, rapidly interconverting conformers,
configurationally unstable stereocenters, ligand exchange in coordination complexes).

For a site s, let X_s be its finite local geometric configuration space.
The existing local permutation groups in stereodescriptors.py define X_s as a
set of permutation cosets/orbits; e.g. a generic tetrahedron has two choices,
while a generic square-planar site has three.

For a **partially specified** assignment sigma at other sites define

    M(sigma, s=x)

to be the stereo graph after setting site s to local state x.

A site s is stereogenic at sigma if there exist x,y in X_s such that

    M(sigma,s=x) is NOT isomorphic to M(sigma,s=y)

**under an isomorphism constrained to map s to itself**, while respecting all
already fixed configurations and chemical labels. Anchoring s is important:
unanchored molecular isomorphism may exchange equivalent site positions.

The fully unrestricted problem asks whether *some* completion of sigma makes
this true. That existential problem explains why an initially equal-color
ligand pair must **not** automatically disqualify a potential site.

## 2. Why the current code is useful

* `StereoMolGraph.atom_stereo` and `StereoMolGraph.bond_stereo` already carry
  geometry-specific descriptors.
* `Stereo.get_isomers()` already enumerates local geometric configurations,
  including non-tetrahedral arrangements.
* `color_refine_mg` gives a fast first-stage topological partition.
* `vf2pp_all_isomorphisms(..., stereo=True)` provides a correctness fallback
  for comparing complete stereo assignments.
* `RDMol2StereoMolGraph` can nominate configurations from RDKit geometry,
  hybridization, and bond tags.

A plain `MolGraph` does **not** in general retain bond order, formal charge,
implicit hydrogens, or coordination geometry. Deriving universal stereogenicity
directly from bare connectivity is not well posed. The *geometry nomination*
step must therefore remain separate from *stereogenicity testing*.

## 3. Proposed architecture

    chemical graph / geometry metadata
                 |
                 v
      geometry candidate nomination
                 |
                 v
      local configuration domains X_s
                 |
                 v
     fast color/ligand partition checks
                 |
                 v
    dependency / influence graph on sites
                 |
                 v
    incremental local-state enumeration
                 |
                 v
    exact completion check + duplicate check
                 |
                 v
     unique, validated stereo graphs

### Phase A: Geometry candidate nomination

Use existing descriptors as the authoritative source of geometry:
`Tetrahedral`, `PlanarBond`, `AtropBond`, `SquarePlanar`,
`TrigonalBipyramidal`, `Octahedral`, and other supplied descriptors.
A specified descriptor is a fixed constraint; an unspecified descriptor
represents an allowed local domain.

For an initial automatic shortcut, consider fully explicit four-coordinate
carbon/silicon atoms as tetrahedral **candidates**. Do **not** automatically
interpret every 3-coordinate atom as tetrahedral: three neighbors can mean
a hidden hydrogen, lone pair, pi bond, or a trigonal geometry.

Future perception should accept an optional chemistry adapter
(e.g. RDKit input with explicit bond orders/charges/implicit H/coordination).
Candidate nomination should be conservative about *rejection*: an unrecognized
shape should be marked unsupported rather than silently discarded.

For every candidate site s, build the local domain X_s once, taking
advantage of the existing permutation-group representation. Cache
canonical local configurations so equivalent permutations are collapsed
*within* the geometry. This local quotient is distinct from pruning
global molecular automorphism orbits.

### Phase B: Fast classification of ligands

Start from atomic labels and 1-dimensional Weisfeiler-Lehman color refinement
of the **connectivity graph without assigning arbitrary stereo parities**.

For a tetrahedron, compare the colors of its four ligand attachments.
For planar bonds and atropisomer bonds, compare each end separately.

* Pairwise distinct attachment colors are a sufficient test that the
  attachments are constitutionally distinguishable. With a supported
  local geometry this provides a fast sufficient condition for local
  stereogenicity (subject to physical feasibility).
* Equal colors do **not** prove equivalent ligands: 1-WL can fail on cyclic
  and highly symmetric graphs.
* Even if attachments are constitutionally equivalent, they may become
  distinguishable after other stereochemical assignments are made.
  Therefore unresolved candidates are retained, never dropped by 1-WL.

The first prototype uses statuses:

* `specified`: a parity/configuration is already fixed;
* `distinct_colors`: all relevant ligand colors differ;
* `unresolved`: detection needs an exact symmetry/assignment analysis.

These statuses are intentionally **not** a complete stereogenicity certificate.

### Phase C: Dependency graph

Let D(s) mean *sites that might influence the stereogenicity of s*.
Edges are directed: t -> s means assignment at t might change whether s
has distinguishable local configurations. An edge is not necessarily
bidirectional.

A **cheap exploratory estimate**, implemented in the prototype, is:

1. Identify color-equal ligand pairs at each unresolved candidate s.
2. Temporarily give the atom(s) defining candidate t a unique seed color.
3. Run color refinement again.
4. If a formerly equal pair at s is split, record t -> s.

This tests **positional influence**, not actual configuration influence.
It can overestimate or underestimate real stereochemical dependencies.
It must **not** be used for correctness-critical pruning.

For a rigorous second version:

* Under fixed partial stereo sigma, compute for each site s the local
  permutations induced by stereo-preserving graph isomorphisms that fix s.
* Search for witness mappings exchanging candidate ligand positions or
  locally equivalent states.
* Record which remote sites and **which local assignments** break each
  witness mapping; simple site-to-site edges are insufficient for
  conjunctive dependencies.
* Represent conditions as a witness-labelled hypergraph or as constraints
  on tuples of sites, not necessarily a simple directed graph.
* For true proof of irrelevance, require that *all* possible state
  completions preserve a suitable local equivalence mapping, or use a
  sound overapproximation of influence.

Illustration:

    s (central candidate)
      /         \
  arm A       arm B
   [t]         [u]

The arms may be constitutionally equivalent. Assigning a different
configuration to t and u can distinguish the arms, activating s.
This is why irreversible per-site filtering at the beginning is unsafe.

### Phase D: Enumeration without automorphism pruning (current scope)

The draft code first constructs X_s for every free candidate and
enumerates the Cartesian product lazily:

    for assignment in X_1 x ... x X_k:
        attach local descriptors to a new StereoMolGraph
        if requested, deduplicate completed graphs by stereo isomorphism
        yield the resulting candidate configuration

No full product is stored in memory, no graph copy is needed before a
complete assignment, and an optional output limit is supported.
The runtime remains proportional to the number of tried assignments,
which may grow exponentially. The unique-output option still stores
already emitted graphs.

**Important:** final duplicate removal is not the McKay automorphism-orbit
pruning proposed for a later version. It does not avoid generating
equivalent *partial* assignments.

For the next incremental DFS version, use the dependency graph to choose
sites with smaller local domains and locally resolved constraints first.
Cache color-refinement results for the fixed skeleton and for changed
stereo neighborhoods. This affects speed but not correctness.

### Phase E: Exact local stereogenicity validation (planned, not yet coded)

To avoid labeling a completely non-stereogenic atom as stereogenic:

1. Given a completed candidate assignment, isolate the site s.
2. Generate alternative local states at s, keeping other states fixed.
3. Test whether two alternatives are distinct using the
   **s-anchored stereo-isomorphism** test.
4. If all alternatives are equivalent, s is not a locally
   stereogenic descriptor in that completed context. Remove or
   mark it accordingly, then re-evaluate coupled sites if required.
5. Only then use exact molecular stereo isomorphism to remove
   duplicate completed structures.

Simply comparing two unanchored molecular graphs is insufficient for
step 3 because it may exchange s with an equivalent site. Applying
step 4 in the presence of mutually coupled/conditionally stereogenic
centers needs a well-defined simultaneous normalization rule, not an
unproven greedy deletion strategy.

An eventual exact implementation should provide a completeness argument
and handle ring closures, para-stereogenicity, local stereochemical
isomerism for geometries with >2 states, and multiple components.

## 4. Prototype API

    from stereomolgraph.experimental.stereoisomer_draft import (
        detect_stereo_elements,
        enumerate_stereo_assignments,
    )

    # g is a StereoMolGraph, preferably from a geometry-aware importer
    analysis = detect_stereo_elements(g, infer_tetrahedral=True)

    for site in analysis.sites:
        print(site.key, site.classification, len(site.states))
        print("Potential influences:", analysis.dependencies[site.key])

    for candidate in enumerate_stereo_assignments(
        g, analysis, unique=True, max_isomers=100
    ):
        print(candidate)

Input descriptors with predefined parity are held fixed. To change one,
explicitly replace it with a descriptor whose parity is None before
candidate generation. The current implementation does not change any
existing public API or existing enumerator.

The output name deliberately says **assignments** rather than
**stereoisomers** until Phase E is implemented and validated.

## 5. Expected computational cost

Let n = number of atoms; m = bonds; k = candidate stereo sites;
d_s = number of local configurations at site s; W(n,m) = cost of
1-WL refinement; I(n,m) = one stereo-isomorphism check.

* Skeleton color refinement: roughly O(W(n,m)).
* Prototype dependency graph: O(k * W(n,m) + k^2 * l^2), where
  l is the small number of ligand slots per site.
* Constructing local domains: small for tetrahedral/E-Z; bounded by
  local coordination geometry (square planar/TBP/octahedral are larger).
* Assignment exploration: O(product_s d_s).
* Each emitted candidate requires a graph copy; unique=True may require
  further expensive exact stereo-isomorphism tests (as in the old code).

This is already efficient in **memory and local-domain generation**, but
**not** an algorithmic solution to exponential-time enumeration.
The intended later acceleration is to avoid unnecessary exact checks
and eventually to incorporate canonical augmentation/orbit pruning.

## 6. Validation milestones

1. **Simple local stereo:** one tetrahedron with 4 different ligands;
   E/Z with distinct substituents on both ends.
2. **Ordinary non-stereogenic sites:** methane/CH2X2 with explicit H;
   generated descriptor should ultimately be absent, not just
   counted once after symmetry deduplication.
3. **Meso and pseudoasymmetry:** 2,3-disubstituted butane,
   inositol, highly symmetric bicyclic structures.
4. **Coupled sites:** two constitutionally identical arms with different
   remote stereo assignments, including rings.
5. **General geometries:** square planar, trigonal bipyramidal,
   octahedral, atropisomer bonds, and missing ligands (None).
6. **Atom relabeling:** site classification, dependency structure and
   enumerated isomer classes must be invariant under renumbering.
7. **Cross-validation:** compare the *number and set* of valid isomers
   against independent chemistry references, RDKit, and enu where
   their scopes overlap. Do not assume agreement proves correctness.
8. **Feasibility:** distinguish graph assignments from geometrically
   realizable and configurationally stable stereoisomers.

## 7. Deliberately omitted

* Full chemistry-aware geometry detection from bare connectivity.
* Certified dependencies and complete stereo-element detection.
* Global automorphism-based branch/orbit pruning.
* Stereochemical assignment normalization for inactive sites.
* Geometry feasibility, kinetics of interconversion, and stereochemical
  constraints of rigid cycles.
* Enantiomer-pair quotienting or CIP priority assignment.
* Multi-step transition-state/reaction-graph stereo enumeration.

The design aims to stabilize the local-state and dependency interfaces
before introducing complexity from canonicalization or McKay pruning.
