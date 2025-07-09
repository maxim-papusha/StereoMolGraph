Stereodescriptors
==================

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
    :special-members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.SquarePlanar
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.TrigonalBipyramidal
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.Octahedral
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.PlanarBond
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.AtropBond
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.AtomStereo
    :members:
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.BondStereo
    :members:
    :inherited-members: