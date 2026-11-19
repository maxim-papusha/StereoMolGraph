# Symmetry

### stereomolgraph.algorithms.symmetry.bond_automorphism_classes(smg: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), mappings: Iterable[Mapping[int, int]] | None = None) → dict[frozenset[int], set[frozenset[int]]]

Bonds are only considered if both of its atoms have other substituents
and therefore can have a internal rotation.

* **Return type:**
  `dict`[`frozenset`[`int`], `set`[`frozenset`[`int`]]]

### stereomolgraph.algorithms.symmetry.symmetry_number(graph: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), atom_labels=None) → int

Calculated from the number of graph isomorphisms which conserve the
stereo information.
symmetry_number = internal_symmetry_number \* rotational_symmetry_number

* **Return type:**
  `int`

### stereomolgraph.algorithms.symmetry.external_symmetry_number(graph: [StereoMolGraph](../graphs/stereomolgraph.md#stereomolgraph.StereoMolGraph), mappings: Iterable[dict[int, int]] | None = None) → int

Calculate the upper bound of the external symmetry number for StereoMolGraph

* **Return type:**
  `int`
