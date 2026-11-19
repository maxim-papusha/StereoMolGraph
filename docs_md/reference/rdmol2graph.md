# RDMol2Graph

### *class* stereomolgraph.rdmol2graph.RDMol2StereoMolGraph(stereo_complete: bool = False, use_atom_map_number: bool = False, lone_pair_stereo: bool = True, resonance: bool = True, \_max_resonance_structures: int = 1000, \_radical_resonance_as_cation: bool = True, \_min_trans_ring_size: int = 7) → None

Convert an RDKit `rdkit.Chem.Mol` to a `StereoMolGraph`.

All aromatic bonds are considered to be cis. Double bonds in rings are
assumed to be cis for rings of size <= 7.

* **Parameters:**
  * **stereo_complete** (`bool`) – If `True`, attempt to infer complete stereo
    parities when possible. Defaults to `False`. (default: `False`)
  * **use_atom_map_number** (`bool`) – Use RDKit atom map numbers as atom identifiers
    instead of RDKit atom indices. Defaults to `False`. (default: `False`)
  * **lone_pair_stereo** (`bool`) – Include stereochemistry that depends on lone
    pairs (if present). Defaults to `True`. (default: `True`)
  * **resonance** (`bool`) – Enumerate resonance structures and merge bond stereo
    information from them. Defaults to `True`. (default: `True`)
