# Isomorphism

### stereomolgraph.algorithms.isomorphism.vf2pp_all_isomorphisms(g1: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph) | [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph) | [CondensedReactionGraph](../graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph) | [StereoCondensedReactionGraph](../graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph), g2: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph) | [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph) | [CondensedReactionGraph](../graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph) | [StereoCondensedReactionGraph](../graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph), atom_labels: None | tuple[np.ndarray[tuple[int], np.dtype[np.int64]], np.ndarray[tuple[int], np.dtype[np.int64]]] = None, stereo: bool = False, stereo_change: bool = False, subgraph: bool = False) → Iterator[dict[AtomId, AtomId]]

Find all isomorphisms between two graphs [[1](../../index.md#id4)].

* **Parameters:**
  * **g1** ([`MolGraph`](../graphs/molgraph.md#stereomolgraph.MolGraph) | [`StereoMolGraph`](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph) | [`CondensedReactionGraph`](../graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph) | [`StereoCondensedReactionGraph`](../graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)) – First graph
  * **g2** ([`MolGraph`](../graphs/molgraph.md#stereomolgraph.MolGraph) | [`StereoMolGraph`](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph) | [`CondensedReactionGraph`](../graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph) | [`StereoCondensedReactionGraph`](../graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)) – Second graph
  * **atom_labels** (`None` | `tuple`[`ndarray`[`tuple`[`int`], `dtype`[`int64`]], `ndarray`[`tuple`[`int`], `dtype`[`int64`]]]) – Optional precomputed graph-ordered atom-label arrays for both
    graphs, if none defaults to color refinement. (default: `None`)
  * **stereo** (`bool`) – Whether to consider stereochemistry in the isomorphism (default: `False`)
  * **stereo_change** (`bool`) – Whether to consider stereochemistry changes in the isomorphism (default: `False`)
  * **subgraph** (`bool`) – Whether to find subgraph isomorphisms instead of graph isomorphisms
    (Currently only limited support and not well tested.) (default: `False`)
* **Return type:**
  `Iterator`[`dict`[`int`, `int`]]
* **Returns:**
  An iterator of all isomorphisms, where each isomorphism is represented as a
  dictionary mapping atom ids of g1 to atom ids of g2.
