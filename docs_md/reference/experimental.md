# Experimental

#### WARNING
The `stereomolgraph.experimental` package contains features whose
interfaces may change without notice. Use at your own risk.

<a id="module-stereomolgraph.experimental"></a>

Experimental functionality for StereoMolGraph.

#### WARNING
This module contains **experimental** features. The interfaces may
change without notice between versions. Use at your own risk.

### stereomolgraph.experimental.generate_fleeting_stereoisomers(graph: [StereoCondensedReactionGraph](graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph), enantiomers: bool = True, atoms: None | Iterable[AtomId] = None, bonds: None | Iterable[Bond] = None) → Iterator[[StereoCondensedReactionGraph](graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)]

Generate all unique fleeting stereoisomers of a
StereoCondensedReactionGraph.

Stereoisomers are constructed by enumerating all combinations of
fleeting stereochemical parities. Only fleeting stereocenters with
undefined parity (`None`) are considered; predefined parities remain
unchanged.

If `enantiomers` is `False`, both enantiomers of each fleeting
stereoisomer are returned. Otherwise, only one representative per
enantiomeric pair is included.

* **Parameters:**
  * **enantiomers** (`bool`) – Whether to remove enantiomers. (default: `True`)
  * **atoms** (`None` | `Iterable`[`int`]) – Optional subset of atoms to consider for stereoisomerism. (default: `None`)
  * **bonds** (`None` | `Iterable`[`frozenset`[`int`]]) – Optional subset of bonds to consider for stereoisomerism. (default: `None`)
* **Yields:**
  Unique fleeting stereoisomers of the input graph.
* **Return type:**
  `Iterator`[[`StereoCondensedReactionGraph`](graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)]

### stereomolgraph.experimental.generate_stereoisomers(graph: [StereoMolGraph](graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), enantiomers: bool = False, atoms: None | Iterable[int] = None, bonds: None | Iterable[frozenset[int]] = None) → Iterator[[StereoMolGraph](graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph)]

Generate all unique stereoisomers of a StereoMolGraph.

Stereoisomers are constructed by enumerating all combinations of
stereochemical parities. Only stereocenters with undefined parity
(`None`) are considered; predefined parities remain unchanged.

If `enantiomers` is `False`, both enantiomers of each stereoisomer
are returned. Otherwise, only one representative per enantiomeric pair
is included.

* **Parameters:**
  * **enantiomers** (`bool`) – Whether to remove enantiomers. (default: `False`)
  * **atoms** (`None` | `Iterable`[`int`]) – Optional subset of atoms to consider for stereoisomerism. (default: `None`)
  * **bonds** (`None` | `Iterable`[`frozenset`[`int`]]) – Optional subset of bonds to consider for stereoisomerism. (default: `None`)
* **Yields:**
  Unique stereoisomers of the input graph.
* **Return type:**
  `Iterator`[[`StereoMolGraph`](graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph)]

### *class* stereomolgraph.experimental.JSONHandler

Serialize and deserialize StereoMolGraph-related graphs to JSON.
