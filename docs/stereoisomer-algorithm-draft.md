# Stereogenicity analysis and stereoisomer enumeration in StereoMolGraph

**Version:** 0.2, research design and initial working prototype  
**Branch:** stereoisomer  
**Implementation:** src/stereomolgraph/experimental/stereoisomer_draft.py

> Scope contract: **The input is always a StereoMolGraph whose existing
> atom and bond stereodescriptors enumerate the complete set of potential
> stereochemical elements.** A descriptor with parity=None has unspecified
> configuration. Detection MUST NOT infer, insert, or remove descriptor sites.
> Enumeration may **replace** an existing undefined descriptor with an
> existing allowed configuration, but may never create a new stereo site.
>
> The current implementation enumerates *stereo-configured graphs*, not
> necessarily physically realizable or fully normalized stereoisomers.
> It does not perform McKay-style automorphism pruning.

## 1. Problem statement and representation

Let the supplied molecular graph be

$$
G=(V,E,\lambda_V,\lambda_E,\mathcal D),
$$

where $V$ is the atom set, $E$ is the bond set, $\lambda_V$ and $\lambda_E$
are the available structural attributes, and

$$
\mathcal D = \{d_s:s\in S\}
$$

is the set of **already present** stereodescriptors. The site key $s$ is an
atom or an unordered pair of bonded atoms, including its kind (atom/bond).
This must be a stable identity across enumeration. No new site is inferred.

Each $d_s$ stores:

- a geometry class (tetrahedral, planar bond, atrop bond, square planar,
  trigonal bipyramidal, octahedral, or a supplied rigid-bond class);
- an ordered tuple of atoms/ligand slots, possibly including the sentinel
  None for a missing ligand;
- a local permutation group $H_s$ that preserves a represented geometry;
- an optional parity, with None meaning **unspecified**, not a wildcard
  in exact equality comparisons.

Some molecules have only undefined descriptors, while some already have
fixed stereo configurations. The algorithm retains the fixed configurations.

**Important modeling boundary:** The constructor of StereoMolGraph or its
RDKit adapter is responsible for nominating geometry-eligible sites and
maintaining any bond-order/chemical context required to justify them.
The algorithm here does not infer geometry from connectivity.

### 1.1 Local state spaces

Let $P_s$ denote the relevant permutations of the ligand positions,
and let $H_s$ denote the geometry-preserving permutation subgroup.
For **distinctly labeled** ligands the local orientation space is
represented schematically by the cosets

$$
X_s \simeq P_s/H_s,\qquad
|X_s| = \frac{|P_s|}{|H_s|}.
\tag{1}
$$

The quotient-side convention depends on whether permutations act on
positions on the left or right; the implementation is defined by
the existing local descriptor equality and permutation conventions.

The familiar generic examples with all ligands distinguishable are:

| Geometry | Local states for distinct ligands |
|---|---:|
| Tetrahedral (4 slots) | $4!/12=2$ |
| Square planar (4 slots) | $4!/8=3$ |
| Trigonal bipyramidal (5 slots) | $5!/6=20$ |
| Octahedral (6 slots) | $6!/24=30$ |
| Planar alkene bond | 2 |
| Axial atropisomer descriptor | 2 (for the supported representation) |

With repeated ligand identities, only a subset of those formal states
might describe distinct local environments. The existing **get_isomers()**
method generates a concrete finite set of descriptors for the given input
site. The implementation uses this set rather than hardcoding binary states.

Define the variable domain of the **existing** descriptor $d_s$ as

$$
D_s =
\begin{cases}
\{d_s\}, & \text{if its parity is already specified},\\[2pt]
\operatorname{get\_isomers}(d_s),
  & \text{if its parity is None}.
\end{cases}
\tag{2}
$$

Some rigid-bond descriptors currently return only themselves, even when
their parity is None. The prototype classifies these as **geometry_only**
rather than incorrectly inventing two configurations.

The initially allowed *assignment space* is

$$
\Omega = \prod_{s\in S}D_s,\qquad
|\Omega|=\prod_{s\in S}|D_s|.
\tag{3}
$$

This is a space of assignments to existing descriptor slots, not yet a
space of rigorously distinct molecular stereoisomers.

### 1.2 Four different meanings of "stereogenic"

We must distinguish:

1. **Candidate:** an input site $s\in S$. The candidates are given.
2. **Clearly distinguishable by constitution:** a site whose different
   ligand attachments receive different suitable graph-refinement colors.
   This provides a fast sufficient distinction for a supported geometry.
3. **Stereogenic at a fixed context:** changing its state while keeping
   all other stereo assignments fixed can change the rooted stereo graph.
4. **Potentially/conditionally stereogenic:** some assignment of *other*
   sites makes it stereogenic. Its status may change with those assignments.

These are neither interchangeable nor simple intrinsic atom properties.
The use of the word "stereoisomer" additionally depends on whether
geometry feasibility and fast interconversion are included in the model.

## 2. Geometry qualification and initial detection

### 2.1 Read-only candidate acquisition

The only candidates come from

~~~
graph.atom_stereo.values()
graph.bond_stereo.values()
~~~

For every supplied candidate, compute $D_s$ once and cache the resulting
local descriptors. In the baseline code, a specified parity is a singleton
domain. No site is added, even for a four-coordinate carbon without a
supplied descriptor. No site is removed on the basis of color equivalence.

A separate candidate-nomination procedure may run *before* this algorithm,
but it is explicitly outside this module and the present specification.

### 2.2 Color refinement

Let $c_0(v)$ be the initial atomic-number color, and for round $r$ define

$$
c_{r+1}(v) =
\operatorname{Compress}\left(
 c_r(v),\
 \operatorname{multiset}\{c_r(u):u\in N(v)\}
\right).
\tag{4}
$$

In the current code, **color_refine_mg** implements a 1-dimensional
Weisfeiler--Lehman-like refinement. The actual implementation uses
integer hashing rather than collision-free symbolic tuples, so numerical
hash collisions are an engineering caveat; color inequality is intended
as a sound structural distinction only conditional on collision-free
refinement and suitably complete structural labels.

Let $L_s=(\ell_1,\ldots,\ell_k)$ be the ligand attachment atoms at site $s$.
A quick sufficient test for fully distinguishable attachments is

$$
c_\infty(\ell_i)\ne c_\infty(\ell_j)
\quad\text{for all }i\ne j.
\tag{5}
$$

For a bond stereo unit, evaluate the ligand pairs at **each end**
independently. The sentinel None is its own color category.
Not all geometry classes have a meaningful pairwise-color test.

**Crucial asymmetry:**

$$
c_\infty(u)\ne c_\infty(v)
\implies
\text{attachments are distinguishable},
$$

but

$$
c_\infty(u)=c_\infty(v)
\centernot\implies
\text{attachments are equivalent}.
\tag{6}
$$

Equal WL colors occur for non-isomorphic rooted environments in sufficiently
complex rings and regular graphs. They can also remain equal before remote
stereochemical configurations distinguish otherwise constitutionally
equivalent branches.

Hence equal colors produce **unresolved**, never **non-stereogenic**.

### 2.3 Prototype site classification

The current **detect_stereo_elements(graph)** returns each input site
with a local domain and one of:

| Status | Meaning |
|---|---|
| specified | Input descriptor has fixed parity; domain is a singleton |
| distinct_colors | Undefined descriptor; quick sufficient ligand-color distinction |
| unresolved | Undefined descriptor; further anchored stereo analysis required |
| geometry_only | Undefined descriptor for which get_isomers() provides only one state |

Only these input sites are returned. The method does not promise that
an unresolved site is stereogenic, nor that it can be discarded.

## 3. Exact local stereogenicity at a fixed assignment

This is the key rigorously defined next layer and has a first
implementation in **locally_stereogenic_at_assignment()**.

Let $\sigma_{\neg s}$ specify every other **variable** site. For each
local state $x\in D_s^{\mathrm{all}}$, form the graph

$$
G_{s,x}\bigl(\sigma_{\neg s}\bigr),
\tag{7}
$$

where $D_s^{\mathrm{all}}$ means all possible local states at $s$
even if its input parity was fixed for the purpose of enumeration.
The other sites remain unchanged.

Let

$$
G_{s,x} \cong_s G_{s,y}
$$

mean an isomorphism that preserves graph attributes, all stereochemical
configurations, and also **fixes site $s$ setwise**:

$$
\phi(s)=s.
\tag{8}
$$

For an atom site, this fixes its central atom. For a bond site, it fixes
the unordered pair of bond endpoints and may exchange them.

A site is locally stereogenic at $\sigma_{\neg s}$ precisely if

$$
\boxed{
A_s(\sigma_{\neg s})=
\begin{cases}
1,&\exists x,y\in D_s^{\mathrm{all}}\colon
G_{s,x}(\sigma_{\neg s})
\not\cong_s
G_{s,y}(\sigma_{\neg s}),\\
0,&\text{otherwise}.
\end{cases}
}
\tag{9}
$$

**Why is the anchoring necessary?** An unrestricted isomorphism could
map the chosen site to a different, equivalent site elsewhere in the
molecule. That answers a global molecular identity question rather
than the local stereogenicity question.

### 3.1 Implementation of the anchored comparison

To anchor an atom site, give that atom a distinguished seed color.
To anchor a bond site, distinguish the two endpoints setwise.
Apply color refinement to both graphs using the **same seed scheme**,
then run the existing stereo-aware VF2++ isomorphism search.

In pseudocode:

~~~python
def locally_stereogenic_at_assignment(G, s):
    require_every_other_variable_site_is_specified(G, s)
    X = all_local_states_of_existing_descriptor(G, s)
    if len(X) < 2:
        return False

    if structurally_distinct_ligands(G, s):
        return True  # fast sufficient condition

    reference = set_existing_descriptor(G.copy(), s, X[0])
    for x in X[1:]:
        alternative = set_existing_descriptor(G.copy(), s, x)
        if not rooted_stereo_isomorphic(reference, alternative, root=s):
            return True
    return False
~~~

This exact test is for **represented stereo graphs**, assuming the stereo
isomorphism routine is correct. It does not resolve unknown physical
geometry, configurational stability, or whether a formally represented
pseudo-descriptor should be removed.

**Current limitations:** The code does not perform this check during
normal enumeration. It is available as a separate diagnostic for
completed assignments (or when only the tested site remains unspecified).
It is not safe to treat parity=None elsewhere as an exact completed
configuration. The old StereoMolGraph equality behavior for
undefined parity must not be used as a wildcard equivalence oracle.

### 3.2 Global potential stereogenicity

A site is potentially stereogenic if it becomes active in at least one
completion of the other local assignments:

$$
P_s =
\max_{\sigma_{\neg s}\in\Omega_{\neg s}}
A_s(\sigma_{\neg s}).
\tag{10}
$$

Calculating this exact existential predicate may itself require
enumerating many remote stereo assignments. The initial prototype
does **not** claim to compute (10); it only classifies candidates
and can test (9) in completed contexts.

In particular, we must **not** reject a site because $A_s=0$ in a
single arbitrary completion. Such rejection loses para-stereogenic
and mutually dependent cases.

## 4. Dependency graph: discovery and semantics

A directed dependency from site $t$ to $s$ means that changing the
configuration at $t$ can alter whether $s$ is locally stereogenic.

The formal edge relation is

$$
t\longrightarrow s
\iff
\exists \tau,\ x,y\in D_t^{\mathrm{all}}:
A_s(\tau,t=x)\ne A_s(\tau,t=y),
\tag{11}
$$

where $\tau$ assigns all relevant sites other than $s,t$.
Equation (11) is an expensive **exact** definition; it is not the
current fast dependency algorithm.

### 4.1 Initial positional influence heuristic

The current implementation follows these steps:

1. Compute unmarked colors $c^{(0)}$ from the input connectivity.
2. For each potential influencing site $t$, give its central atom
   (or both bond endpoints) a special marker in the seed labels.
3. Recompute colors $c^{(t)}$.
4. For each unresolved site $s$, find any two attachment atoms that
   had equal $c^{(0)}$ but different $c^{(t)}$ colors.
5. If such a pair exists, report $t\to s$ in the heuristic dependency map.

The corresponding predicate is

$$
D^{\mathrm{mark}}_{t,s} =
\mathbf 1\left[
\exists u,v\in L_s:
c^{(0)}(u)=c^{(0)}(v)
\land
c^{(t)}(u)\ne c^{(t)}(v)
\right].
\tag{12}
$$

This is a **positional influence detector**, not a proof that
the actual configuration of $t$ influences $s$. It can have false
positives or false negatives, so it must never be used to prune.

A site can depend on several remote sites simultaneously.
For two constitutionally equivalent branches ending in sites $t$ and $u$,
the status at a central site $s$ may be governed by a joint condition
$f(\sigma_t,\sigma_u)$ rather than independent edges. A future
certified dependency representation may therefore need **hyperedges**
or witness-labelled logical constraints.

The current graph is useful for visualization, diagnostics,
and eventual *heuristic ordering* only.

## 5. Initial working enumeration: no automorphism pruning

Enumerate all assignments in (3) lazily, keeping already specified
descriptors fixed:

~~~python
analysis = detect_stereo_elements(G)

for assignment in product(*(site.states for site in analysis.sites)):
    candidate = copy(G)
    for (site, state) in zip(analysis.sites, assignment):
        assert site.key in G.stereo
        replace_descriptor_at_existing_site(candidate, site.key, state)
    if unique and candidate_is_isomorphic_to_previously_seen(candidate):
        continue
    yield candidate
~~~

The code uses the original **StereoMolGraph.set_atom_stereo()** and
**set_bond_stereo()** to replace states **at the same keys**.
It checks that the analysis-site set equals the graph descriptor-site set.
There is no inference or insertion of new sites.

With **unique=True**, completed stereo graphs are compared through the
existing graph hash and stereo-aware isomorphism equality. This is
**leaf-level deduplication**, not automorphism pruning of the search tree.

Its candidate-generation complexity is bounded by

$$
N_\mathrm{attempted}=
\prod_{s\in S}|D_s|
\tag{13}
$$

before global deduplication. The explicit result count may itself be
exponential in the number of genuinely independent stereo units.

Memory usage excluding the stored deduplication set is roughly

$$
O(|V|+|E|+\sum_{s\in S}|D_s|)
\tag{14}
$$

per produced graph, assuming the local state domains have been built.

For comparison, if $N$ candidates are produced and the deduplication
set retains up to $N$ graphs, total memory can grow as $O(N|G|)$.
A future streaming canonical key can reduce the constant factors.

### 5.1 Why this is a correct baseline for represented assignments

For fixed input $S$ and finite local domains $D_s$, the Cartesian
product visits each assignment in $\Omega$ exactly once as a tuple.
Each assignment overwrites *only* the existing descriptors. Therefore:

- no permitted local combination is skipped;
- no unknown descriptor site is added;
- the input graph is not mutated;
- each completion can be checked against an exact stereo graph-isomorphism
  oracle for global duplicates.

This proves exhaustiveness of **represented local assignments**,
not completeness with respect to all physically stable stereoisomers.

When unique=True, correctness of the distinct-output set additionally
depends on the correctness and hash consistency of graph equality,
including behavior for undefined/geometry-only descriptors.

### 5.2 Why global site activation cannot be handled greedily

For a completed assignment, a site can appear locally non-stereogenic
but may become stereogenic after changing another site's configuration.
Simply permanently discarding such a site while exploring partial
assignments would make the enumeration incomplete.

Also, two redundant stereo descriptors can interact: deleting one may
change the apparent status of another. A canonical normalization of
inactive descriptors should be defined and proved, not introduced as
an unverified per-site cleanup rule.

For now the output function is deliberately named
**enumerate_stereo_assignments**, not **enumerate_stereoisomers**.

## 6. Incremental DFS: useful next optimization without group pruning

Instead of the Cartesian product, enumerate states using DFS with
reversible assignments. Let $\sigma_j$ be the partial configuration of
the first $j$ selected sites:

$$
\sigma_0=\varnothing,\qquad
\sigma_{j+1}=\sigma_j\cup\{s_j\mapsto x\}.
\tag{15}
$$

Branching order can favor smaller domains $|D_s|$, fixed sites, and
sites with high dependency-graph connectivity. These are only heuristic
ordering choices and cannot change which leaves are visited.

Fast incremental data can include:

- a cached base color partition;
- local stereo descriptor domain lookups;
- stereo-aware color signatures for the already specified sites;
- a list of unresolved sites affected by the latest configuration.

**Do not** prune a branch merely because an unresolved site has equal
colors or appears inactive before its dependencies are assigned.
That is an unsound shortcut for conditional stereogenicity.

The initial product iterator already has low *candidate-list* memory
requirements; DFS would mainly enable future incremental updates.

## 7. Optional future McKay-like symmetry pruning

This is **not needed for the working version** and is not currently
implemented. Nevertheless it fits naturally with the input descriptor
representation if the group action is carefully defined.

Let

$$
A = \operatorname{Aut}(G_0;\text{attributes and fixed stereo}),
\tag{16}
$$

where $G_0$ is the structural graph carrying all candidate
geometry types but *without arbitrarily specified variable parities*.
Every $g\in A$ maps stereo sites onto compatible stereo sites and
induces a **transport action** on their local configurations.

Writing $g(s)$ for the target site and $g_\ast$ for the transported
local state, the action on full assignments is

$$
(g\cdot\sigma)_{g(s)} = g_\ast(\sigma_s).
\tag{17}
$$

The desired molecular configuration classes are the orbits

$$
\Omega / A,
\tag{18}
$$

for a correctly defined group action and fixed-site constraints.

### 7.1 Safe, small first symmetry optimization

Suppose a DFS prefix $\sigma_j$ is fixed and we are about to assign
site $s_j$. Define the stabilizer of both the prefix and site:

$$
H_j =
\operatorname{Stab}_A(\sigma_j)
\cap \operatorname{Stab}_A(s_j).
\tag{19}
$$

Elements of $H_j$ preserve the already assigned prefix and fix the
current site setwise. They act on its local state set $D_{s_j}$.
It is safe to branch on just one representative of each orbit

$$
D_{s_j}/H_j
\tag{20}
$$

**provided** the complete candidate domains and fixed constraints are
preserved under the action and all remaining assignments are still
enumerated consistently.

This reduces redundant sibling branches but does not resolve all
equivalences arising from permutations of *different unassigned sites*.
It is nonetheless a useful step beyond leaf-level deduplication.

The difficult requirement is a correct transporter
$g_\ast$ for **all** existing descriptor classes:
tetrahedral, planar bonds, atropisomeric, square-planar, TBP,
octahedral, and rigid-bond representations.
It must handle site exchange, bond-endpoint exchange, repeated ligands,
parity inversion and ordered positions without guessing.
One should test the action law

$$
(gh)_\ast(x)=g_\ast(h_\ast(x))
\tag{21}
$$

and verify that every transformed domain remains valid.

### 7.2 Full canonical augmentation

A stronger McKay-like method would use an equivariant
**canonical parent** operation on partially specified stereo graphs.

Given a child partial assignment $C$, compute a canonical
representative (or canonical labelling) and a canonical parent
$p_{\rm can}(C)$. Accept the augmentation from parent $P$ only if it
is equivalent to this canonical parent:

$$
\operatorname{Accept}(P\to C)
\iff
P \cong p_{\rm can}(C).
\tag{22}
$$

A well-defined canonical-parent operation and an appropriate
equivalence criterion are essential. Merely comparing an arbitrary
canonical hash after a choice is not a proof that entire rejected
subtrees are redundant.

The existing **canon** branch offers refinement and backtracking
components for this future work, but it does not currently provide
the complete group generators, stabilizer chain and stereodescriptor
transport required for a safe implementation.

**Recommended gate:** add symmetry pruning only after exhaustive small
cases agree exactly with the unpruned baseline, including all generated
stereo graphs and not just their counts.

## 8. Computational scaling

Let $n=|V|$, $m=|E|$, $k=|S|$, $W(n,m)$ the cost of 1-WL refinement,
and $I(n,m)$ an exact stereo isomorphism cost.

| Stage | Current cost, schematic |
|---|---|
| Read input descriptor slots | $O(k)$ |
| Generate local domains | $O(\sum_s |D_s|)$ plus local permutation work |
| Initial color refinement | $W(n,m)$ |
| Heuristic influence graph | $O(kW(n,m)+k^2L^2)$, $L$ max local slots |
| Generate candidate assignments | $\Theta(\prod_s |D_s|)$ |
| Unique output | Costs of hashing and exact stereo-isomorphism tests |
| Exact local test at one completed site | At most $(|D_s^{\rm all}|-1)$ rooted tests |
| McKay pruning | Not implemented |

No algorithm can explicitly output $2^k$ genuinely distinct binary
stereoisomers in sub-exponential **total output time**. A useful speed
claim must be phrased relative to the number of inequivalent outputs
and the amount of redundant search avoided.

## 9. API and current code

~~~python
from stereomolgraph.experimental.stereoisomer_draft import (
    detect_stereo_elements,
    enumerate_stereo_assignments,
    locally_stereogenic_at_assignment,
)

# G is a StereoMolGraph. ALL candidate descriptors are already stored.
analysis = detect_stereo_elements(G)

for site in analysis.sites:
    print(site.key, site.classification, len(site.states))
    print("Heuristic influences:", analysis.dependencies[site.key])

for completed in enumerate_stereo_assignments(
    G, analysis, unique=True, max_isomers=100
):
    # Site-anchored exact local test, only after all other variable sites
    # are specified (otherwise it raises ValueError).
    for site in analysis.sites:
        active = locally_stereogenic_at_assignment(completed, site.key)
        print(site.key, active)
~~~

This API has no geometry inference option. It does not change the
pre-existing **generate_stereoisomers** public function.

## 10. Verification plan and open questions

### Minimum working regression suite

1. No input descriptors -> no descriptor inserted, one unchanged graph.
2. One supplied undefined tetrahedron with four different ligands:
   exactly two local assignments; exact local test returns true.
3. One supplied undefined tetrahedron with repeated ligands:
   two *formal* assignments; exact local test returns false.
4. Two or more supplied sites: all combinations enumerated.
5. Undefined planar bond and supplied non-tetrahedral geometry:
   configurations generated from get_isomers() without inventing sites.
6. Partially specified input: fixed parity never overwritten.
7. Relabeling atom IDs: same candidate count and same global isomer
   classes, modulo relabeling.
8. Chiral rings and pseudoasymmetric examples: compare against RDKit,
   original chemical nomenclature, and independently enumerated cases.
9. Cyclic and mutual dependencies: verify that no candidate is
   eliminated solely by equal initial colors.
10. Full correctness test: compare any future pruning algorithm to
    the current unpruned baseline for every case small enough to exhaust.

### Known limitations of this version

- Input descriptors may describe geometry that is not chemically or
  physically realizable; this stage trusts the input.
- A *complete* exact classification of globally potentially stereogenic
  sites requires a search over remote assignments (Eq. 10).
- The influence graph is heuristic, not a certificate for pruning.
- Leaf-level graph-isomorphism deduplication may retain configurations
  with inactive stereodescriptors rather than normalizing them.
- No physical ring-closure constraints, conformational barriers, or
  fast-interconversion screening are applied.
- No McKay canonical augmentation or automorphism-group pruning is in
  the baseline; an exact descriptor-transport action is a prerequisite.

The first research milestone is a validated baseline that never silently
skips combinations of **pre-existing** stereo configurations.
The second is exact stereogenicity classification and normalization,
and the third is certified group-theoretical pruning.
