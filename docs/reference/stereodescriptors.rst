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


.. autoclass:: stereomolgraph.stereodescriptors.StereoProtocol
    :no-index-entry:
    :members:

.. autoclass:: stereomolgraph.stereodescriptors.AtomStereo
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:

.. autoclass:: stereomolgraph.stereodescriptors.BondStereo
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:

.. autoclass:: stereomolgraph.stereodescriptors.Tetrahedral
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.SquarePlanar
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.TrigonalBipyramidal
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.Octahedral
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.PlanarBond
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.AtropBond
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 


.. autoclass:: stereomolgraph.stereodescriptors._StereoMixin
    :members:
    :member-order: alphabetical

.. autoclass:: stereomolgraph.stereodescriptors._AchiralStereoMixin
    :members:
    :member-order: alphabetical

.. autoclass:: stereomolgraph.stereodescriptors._ChiralStereoMixin
    :members:
    :member-order: alphabetical