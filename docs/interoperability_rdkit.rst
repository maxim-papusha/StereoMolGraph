RDKit Interoperability
=======================

StereoMolGraph is designed to work seamlessly with `RDKit <https://www.rdkit.org/>`_, a popular cheminformatics toolkit.  
`rdmol2graph` and `graph2rdmol` functions allow for easy conversion between StereoMolGraph's graph representations and RDKit's `Mol` objects.


Included Properties
--------------------
- Atomic connectivity (atom indices and bonds).  
- Local stereochemistry (e.g., tetrahedral chirality).
  - only **Explicit Hydrogens**: StereoMolGraph uses only explicit hydrogens! : Hydrogen atoms are always treated as explicit and must be added if required for RDKit interoperability.  


Not included Properties
-------------------------
StereoMolGraph does **not** store:  
- Bond orders (single/double/triple).  
- Formal atomic charges.  
- Aromaticity or global molecular properties.  

RDKit can compute these missing properties from the connectivity and stereochemistry provided by StereoMolGraph when needed.  