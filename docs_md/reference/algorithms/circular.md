# Circular Algorithm

### stereomolgraph.algorithms.circular.numpy_int_multiset_hash(arr: np.ndarray[tuple[int, ...], np.dtype[np.int64]], out: None | np.ndarray[tuple[Literal[1], ...], np.dtype[np.int64]] = None) → np.ndarray

Hash function for a multiset (order-independent with duplicates) of integers.
Works by sorting the elements and then applying the tuple hashing function.

* **Return type:**
  `ndarray`

### stereomolgraph.algorithms.circular.label_hash(mg: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph), atom_labels: Collection[str] = ('atom_type',)) → ndarray[tuple[int], dtype[int64]]

Generates a hash for each atom based on choosen attributes.

* **Parameters:**
  * **mg** ([`MolGraph`](../graphs/molgraph.md#stereomolgraph.MolGraph)) – MolGraph object containing the atoms.
  * **atom_labels** (`Collection`[`str`]) – Iterable of attribute names to use for hashing. (default: `('atom_type',)`)
* **Return type:**
  `ndarray`[`tuple`[`int`], `dtype`[`int64`]]

### stereomolgraph.algorithms.circular.circular_generator(mg: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph), atom_labels: None | ndarray[tuple[int], dtype[int64]] = None) → Iterator[ndarray[tuple[int], dtype[int64]]]

Color refinement algorithm for MolGraph.

This algorithm refines the atom coloring based on their connectivity.
Identical to the Weisfeiler-Lehman (1-WL) algorithm.

* **Parameters:**
  * **mg** ([`MolGraph`](../graphs/molgraph.md#stereomolgraph.MolGraph)) – MolGraph object containing the atoms and their connectivity.
  * **max_iter** – Maximum number of iterations for refinement.
    Default is None, which means it will run until convergence.
* **Return type:**
  `Iterator`[`ndarray`[`tuple`[`int`], `dtype`[`int64`]]]

### stereomolgraph.algorithms.circular.circular_stereo_generator(smg: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), atom_labels: None | ndarray[tuple[int], dtype[int64]] = None) → Iterator[ndarray[tuple[int], dtype[int64]]]

Generator of the Circular Stereo Algorithm for StereoMolGraph [[2](../../index.md#id5)].

* **Return type:**
  `Iterator`[`ndarray`[`tuple`[`int`], `dtype`[`int64`]]]

### stereomolgraph.algorithms.circular.color_refine_hash_mg(graph: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph)) → int

Color-refined hash for plain MolGraph objects.

* **Return type:**
  `int`

### stereomolgraph.algorithms.circular.color_refine_hash_smg(graph: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph)) → int

Color-refined hash for StereoMolGraph objects.

Drops the extra sentinel slot the stereo generator appends.

* **Return type:**
  `int`

### stereomolgraph.algorithms.circular.color_refine_hash_crg(graph: [CondensedReactionGraph](../graphs/condensedreactiongraph.md#stereomolgraph.CondensedReactionGraph)) → int

Color-refined hash for CondensedReactionGraph objects.

* **Return type:**
  `int`

### stereomolgraph.algorithms.circular.color_refine_hash_scrg(graph: [StereoCondensedReactionGraph](../graphs/stereocondensedreactiongraph.md#stereomolgraph.StereoCondensedReactionGraph)) → int

Color-refined hash for StereoCondensedReactionGraph objects.

* **Return type:**
  `int`

### stereomolgraph.algorithms.circular.circular_fingerprint(graph: [MolGraph](../graphs/molgraph.md#stereomolgraph.MolGraph), radius: int = 3, n_bits: int = 2048, count: bool = True, accumulate: bool = True, include_hydrogens: bool = False) → ndarray

Build a circular fingerprint for a molecular graph.

If `accumulate` is `True`, identifiers from every radius up to the
requested radius are included. Otherwise, only the last radius is used.
If `include_hydrogens` is `False`, atom environments centered on
hydrogens are excluded.

* **Parameters:**
  * **graph** ([`MolGraph`](../graphs/molgraph.md#stereomolgraph.MolGraph)) – Molecular graph to fingerprint.
  * **radius** (`int`) – Maximum refinement radius to include. (default: `3`)
  * **n_bits** (`int`) – Length of the folded fingerprint. If `0` or `None`,
    return the unique integer identifiers instead of a folded bit/count vector. (default: `2048`)
  * **count** (`bool`) – If `True`, all unique values are used once. (default: `True`)
  * **accumulate** (`bool`) – Whether to include identifiers from all radii up to
    `radius` instead of only the final radius. (default: `True`)
  * **include_hydrogens** (`bool`) – Whether to include hydrogen-centered
    environments. (default: `False`)
* **Return type:**
  `ndarray`

### stereomolgraph.algorithms.circular.circular_stereo_fingerprint(graph: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), radius: int = 3, n_bits: int = 2048, count: bool = True, accumulate: bool = True, include_hydrogens: bool = False) → ndarray

Build a circular stereo fingerprint for a molecular graph.

If `accumulate` is `True`, identifiers from every radius up to the
requested radius are included. Otherwise, only the last radius is used.
If `include_hydrogens` is `False`, atom environments centered on
hydrogens are excluded.

* **Parameters:**
  * **graph** ([`StereoMolGraph`](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph)) – Molecular graph to fingerprint.
  * **radius** (`int`) – Maximum refinement radius to include. (default: `3`)
  * **n_bits** (`int`) – Length of the folded fingerprint. If `0` or `None`,
    return the unique integer identifiers instead of a folded bit/count vector. (default: `2048`)
  * **count** (`bool`) – If `True`, all unique values are used once. (default: `True`)
  * **accumulate** (`bool`) – Whether to include identifiers from all radii up to
    `radius` instead of only the final radius. (default: `True`)
  * **include_hydrogens** (`bool`) – Whether to include hydrogen-centered
    environments. (default: `False`)
* **Return type:**
  `ndarray`
