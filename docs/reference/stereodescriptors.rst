Stereodescriptors
=================

Stereodescriptors are standardized notations that precisely define the 3D spatial arrangement of atoms around a molecular center. These descriptors combine three key aspects of molecular geometry:

- **Hybridization state**
- **Stereochemical configuration**
- **Relative positioning** of substituents

They serve as fundamental building blocks for representing molecular stereochemical perception in algorithms.

Core Principles
---------------

**Order Sensitivity**  
Stereodescriptor assignments depend critically on the input order of ligand atoms.

**Permutation Invariance**  
Within their defined symmetry group, stereodescriptors are invariant to permitted atom permutations. For example:  


Parity System Explained
-----------------------

The parity value encodes the spatial configuration relative to the input atom order:

+----------+-----------------+------------------------------------------------------------+
| Parity   | Type            | Meaning & Examples                                         |
+==========+=================+============================================================+
| ``1``    | Chiral          | One enantiomer configuration|
|          |                 |                    |
+----------+-----------------+------------------------------------------------------------+
| ``-1``   | Chiral          | |
+----------+-----------------+------------------------------------------------------------+
| ``0``    | Achiral         | Non-chiral center:                                        |
|          |                 |               |
|          |                 |                                    |
+----------+-----------------+------------------------------------------------------------+
| ``None`` | Undefined/Any   | Stereochemical wildcard:                                  |
|          |                 | - Placeholder for unknown stereochemistry                 |
|          |                 | - Permits arbitrary configurations                       |
|          |                 | - Disables stereochemical validation                      |
+----------+-----------------+------------------------------------------------------------+

Practical Implementation
------------------------

**Chiral Centers**  
Tetrahedral atoms and other chiral elements use ``1``/``-1`` parity to distinguish enantiomers.   

**Achiral Centers**  
Use parity ``0``   

**Special Case: ``None``**  
The ``None`` parity is used for:  
1. Represents unknown or unspecified stereochemistry  
2. Allows non-enantiomeric rearrangements (e.g., square planar ligand exchange)  
3. Functions as a stereochemical wildcard in substructure searching    
        

When using stereodescriptors, remember they encode *local* configuration. Global molecular chirality emerges from combinations of local descriptors and molecular topology.


.. autoclass:: stereomolgraph.stereodescriptors.Tetrahedral
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.SquarePlanar
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.TrigonalBipyramidal
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.Octahedral
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.PlanarBond
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.AtropBond
    :members: