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


.. autoclass:: stereomolgraph.stereodescriptors.Stereo
    :no-index-entry:
    :members:
    :member-order: bysource


.. autoclass:: stereomolgraph.stereodescriptors.AtomStereo
    :show-inheritance:
    :no-index-entry:
    :members:
    :member-order: alphabetical

.. autoclass:: stereomolgraph.stereodescriptors.BondStereo
    :show-inheritance:
    :no-index-entry:
    :members:
    :member-order: alphabetical


.. autoclass:: stereomolgraph.stereodescriptors.Tetrahedral
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance: 

.. autoclass:: stereomolgraph.stereodescriptors.Allene
    :no-index-entry:
    :members:
    :member-order: alphabetical
    :inherited-members:
    :show-inheritance:

**Allene atom stereochemistry**

An allene descriptor is stored at the central atom of a C=C=C unit.
Its two terminal atoms are one bond away; the four substituents whose
order defines axial stereochemistry are two bonds away. Consequently,
the circular stereo algorithm first distinguishes its two orientations
at **aggregation radius 2**, not radius 1.

RDKit objects can represent this with the central atom's
``CHI_ALLENE`` tag and ``_chiralPermutation`` (1 or 2 for the
two orientations, or 0 for unspecified stereochemistry). This
conversion conserves atom stereo in RDKit molecule objects, but
standard RDKit SMILES export does not reliably retain allene chirality.
The descriptor's parity is relative and is not an R/S CIP assignment.


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
