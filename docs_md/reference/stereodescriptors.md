# Stereodescriptors

Stereodescriptors are standardized notations that precisely define the 3D spatial arrangement of atoms around a molecular center. These descriptors combine three key aspects of molecular geometry:

- **Hybridization state**
- **Stereochemical configuration**
- **Relative positioning** of substituents

They serve as fundamental building blocks for representing molecular stereochemical perception in algorithms.

**Order Sensitivity**
Stereodescriptor assignments depend critically on the input order of ligand atoms.

**Permutation Invariance**
Within their defined symmetry group, stereodescriptors are invariant to permitted atom permutations.

**Chiral Centers**
Tetrahedral atoms and other chiral elements use `1`/`-1` parity to distinguish enantiomers.

**Achiral Centers**
Use parity `0`

**Special Case: \`\`None\`\`**
The `None` parity is used for:
1. Represents unknown or unspecified stereochemistry
2. Allows non-enantiomeric rearrangements (e.g., square planar ligand exchange)
3. Functions as a stereochemical wildcard in substructure searching

When using stereodescriptors, remember they encode *local* configuration. Global molecular chirality emerges from combinations of local descriptors and molecular topology.

### *class* stereomolgraph.stereodescriptors.Stereo(atoms: A, parity: P = None)

Protocol to represent the orientation of a group of atoms in space.
This is used to represent local stereochemistry and simultaneously the
hybridization of atoms.

#### atoms *: `TypeVar`(`A`, bound= `tuple`[`None` | `int`, `...`], covariant=True)*

Atoms are a order dependent tuple of integers.

#### parity *: `TypeVar`(`P`, bound= `Optional`[`Literal`[`1`, `0`, `-1`]], covariant=True)*

parity is a number that defines the orientation of the atoms. If None,
the relative orientation of the atoms is not defined.
If 0 the orientation is defined and part of a achiral stereochemistry.
If 1 or -1 the orientation is defined and part of a chiral stereochemistry.

#### *property* PERMUTATION_GROUP *: Iterable[A]*

Defines all allowed permutations defined by the symmetry group under
which the stereochemistry is invariant.

#### invert() → Self

Inverts the stereo. If the stereo is achiral, it returns itself.

* **Return type:**
  `Self`

#### get_isomers() → Set[Self]

Returns all stereoisomers of the stereochemistry. Not just the
inverted ones, but all possible stereoisomers.

* **Return type:**
  `Set`[`Self`]

### *class* stereomolgraph.stereodescriptors.AtomStereo(atoms: A, parity: P = None)

Bases: [`Stereo`](#stereomolgraph.stereodescriptors.Stereo)[`A`, `P`], `Protocol`, `Generic`[`A`, `P`]

### *class* stereomolgraph.stereodescriptors.BondStereo(atoms: A, parity: P = None)

Bases: [`Stereo`](#stereomolgraph.stereodescriptors.Stereo)[`A`, `P`], `Protocol`, `Generic`[`A`, `P`]

### *class* stereomolgraph.stereodescriptors.Tetrahedral(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[1, -1]]

Represents all possible configurations of atoms for a Tetrahedral
Stereochemistry:

```default
parity = 1      parity = -1
    4                4
    |                |
    0                0
 /  ¦  \          /  ¦  \
2   1   3        3   1   2
```

Atoms of the tetrahedral stereochemistry are ordered in a way that when the
first atom is rotated to the back, the other atoms in order are rotated in
the direction defined by the stereo.

* **Variables:**
  * **atoms** – Atoms of the stereochemistry
  * **parity** – Stereochemistry
  * **PERMUTATION_GROUP** – Permutations allowed by the stereochemistry

### *class* stereomolgraph.stereodescriptors.SquarePlanar(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[0]]

Represents all possible configurations of atoms for a
SquarePlanar Stereochemistry:

```default
1     4
 \   /
   0
 /   \
2     3
```

Atoms of the Square Planar stereochemistry are ordered in a way that

* **Variables:**
  * **atoms** – Atoms of the stereochemistry
  * **parity** – Stereochemistry

### *class* stereomolgraph.stereodescriptors.TrigonalBipyramidal(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[1, -1]]

Represents all possible configurations of atoms for a
TrigonalBipyramidal Stereochemistry:

```default
parity = 1             parity = -1
 3   1                     1   3
  ◁  ¦                    ¦  ▷
     0  — 5           5 —  0
  ◀  ¦                    ¦  ▶
 4   2                     2   4
```

Atoms of the trigonal bipyramidal stereochemistry are ordered in a way that
when the first two atoms are the top and bottom of the bipyramid. The last
three equatorial atoms are ordered in a way that when the first atom is
rotated to the back, the other atoms in order are rotated in the direction
defined by the stereo.

* **Variables:**
  * **atoms** – Atoms of the stereochemistry
  * **parity** – Stereochemistry

### *class* stereomolgraph.stereodescriptors.Octahedral(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[1, -1]]

Represents all possible configurations of atoms for a Octahedral
Stereochemistry:

```default
parity = 1             parity = -1
 3  1   6                3  2  6
  ◁ ¦ /                  ◁ ¦ /
    0                       0
  / ¦ ▶                  / ¦  ▶
 4  2  5                4   1  5
```

### *class* stereomolgraph.stereodescriptors.PlanarBond(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`None` | `int`, `None` | `int`, `int`, `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[0]]

Represents all possible configurations of atoms for a
Planar Structure and should be used for aromatic and double bonds:

```default
0        4
 \      /
  2 == 3
 /      \
1        5
```

All atoms of the double bond are in one plane. Atoms 2 and 3 are the center
Atoms 0 and 1 are bonded to 2 and atoms 4 and 5 are bonded to 3.
The stereochemistry is defined by the relative orientation
of the atoms 0, 1, 4 and 5.

* **Variables:**
  * **atoms** – Atoms of the stereochemistry
  * **parity** – Stereochemistry
  * **PERMUTATION_GROUP** – Permutations allowed by the stereochemistry

### *class* stereomolgraph.stereodescriptors.AtropBond(atoms: A, parity: P = None)

Bases: `_StereoMixin`[`tuple`[`None` | `int`, `None` | `int`, `int`, `int`, `None` | `int`, `None` | `int`], `None` | `Literal`[1, -1]]

Represents all possible configurations of atoms for a
Atropostereoisomer bond:

```default
parity = 1          parity = -1
1       5           1        5
 \     /            ◀      /
  2 - 3               2 - 3
◀      \            /      \
0        4         0         4
```
