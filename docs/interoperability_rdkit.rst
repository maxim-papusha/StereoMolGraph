RDKit Interoperability
=======================

StereoMolGraph is designed to work seamlessly with `RDKit <https://www.rdkit.org/>`_, a popular cheminformatics toolkit.  

Key Features
-------------
- **Bidirectional Conversion**: All ``Graph`` objects in StereoMolGraph can be converted to RDKit ``Mol`` objects (and vice versa), enabling easy integration into chemoinformatic workflows.  
- **Focused Scope**: StereoMolGraph **only** stores:  
  - Atomic connectivity (atom indices and bonds).  
  - Local stereochemistry (e.g., tetrahedral chirality).
- **Explicit Hydrogens**: Hydrogen atoms are always treated as explicit and must be added if required for RDKit interoperability.  
- Additional properties can be stored in StereoMolGraph objects as attributes atoms or bonds, allowing for custom data to be associated with the graph structure.


Not included Properties
-------------------------
StereoMolGraph does **not** store:  
- Bond orders (single/double/triple).  
- Formal atomic charges.  
- Aromaticity or global molecular properties.  

RDKit can compute these missing properties from the connectivity and stereochemistry provided by StereoMolGraph when needed.  