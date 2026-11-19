# Coords

### *class* stereomolgraph.coords.GeometryProtocol(\*args, \*\*kwargs)

Protocol for geometry-like objects used throughout StereoMolGraph.

Implementations must provide atom types and Cartesian coordinates in the
same atom order.

* **Variables:**
  * **atom_types** – Atom types for each atom in the geometry.
  * **coords** – Cartesian coordinates with shape `(n_atoms, 3)` in
    Angstrom.

### stereomolgraph.coords.are_planar(points: ndarray[tuple[N, Literal[3]], NP_FLOAT], threshold: float = 0.35) → bool

Checks if all atoms are in one plane

Checks if the all atoms are planar within a given threshold.
The threshold is the maximal distance of an atom from the plane of
three other atoms.

* **Parameters:**
  * **points** (*np.ndarray*) – coordinates of the atoms
  * **threshold** (*float*) – maximal distance of atom from plane [Angstrom] (default: `0.35`)
* **Returns:**
  True if all atoms are planar
* **Return type:**
  bool

### stereomolgraph.coords.handedness(coords: ndarray[tuple[Literal[4], Literal[3]], NP_FLOAT]) → ndarray[tuple[Literal[1]], dtype[int8]]

Calculates the orientation of the atom 4 from the plane defined
by the first three atoms from their coordinates.

* **Return type:**
  `ndarray`[`tuple`[`Literal`[`1`]], `dtype`[`int8`]]

### *class* stereomolgraph.coords.Geometry(atom_types: Sequence[int | str] = (), coords: ndarray[tuple[int, int], NP_FLOAT] = array([], shape=(0, 3), dtype=float64))

Represents a molecular geometry, i.e. the coordinates and atom types.

* **Parameters:**
  * **atom_types** (`Sequence`[`int` | `str`]) – tuple of Element objects (default: `()`)
  * **coords** (`ndarray`[`tuple`[`int`, `int`], `TypeVar`(`NP_FLOAT`, bound= `dtype`[`floating`], contravariant=True)]) – Cartesian coordinates of atoms in Angstrom (default: `array([], shape=(0, 3), dtype=float64)`)

#### *classmethod* from_xyz_file(path: PathLike[str] | str) → [Geometry](#stereomolgraph.coords.Geometry)

Create a Geometry from an XYZ file.

* **Return type:**
  [`Geometry`](#stereomolgraph.coords.Geometry)

#### *classmethod* from_xyz(xyz_string: str) → [Geometry](#stereomolgraph.coords.Geometry)

Create a Geometry from an XYZ-format string.

This mirrors [`from_xyz_file()`](#stereomolgraph.coords.Geometry.from_xyz_file) but reads from a string instead of a
file path.

* **Return type:**
  [`Geometry`](#stereomolgraph.coords.Geometry)

#### *classmethod* from_xyz_file_multi(path: PathLike[str] | str) → tuple[[Geometry](#stereomolgraph.coords.Geometry), ...]

Create multiple Geometries from an XYZ file containing multiple structures.

* **Parameters:**
  **path** (`PathLike`[`str`] | `str`) – Path to XYZ file containing one or more structures
* **Return type:**
  `tuple`[[`Geometry`](#stereomolgraph.coords.Geometry), `...`]
* **Returns:**
  Tuple of Geometry objects

#### *classmethod* from_xyz_multi(xyz_string: str) → tuple[[Geometry](#stereomolgraph.coords.Geometry), ...]

Create multiple Geometries from an XYZ-format string containing multiple structures.

* **Parameters:**
  **xyz_string** (`str`) – XYZ-format string containing one or more structures
* **Return type:**
  `tuple`[[`Geometry`](#stereomolgraph.coords.Geometry), `...`]
* **Returns:**
  Tuple of Geometry objects

#### xyz_str(comment: None | str = None) → str

returns the xyz representation of this geometry as a string

* **Parameters:**
  **comment** (`None` | `str`) – comment for 2nd line of xyz file (default: `None`)
* **Return type:**
  `str`
* **Returns:**
  xyz representation of this geometry
