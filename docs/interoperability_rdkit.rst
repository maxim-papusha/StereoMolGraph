RDKit Interoperability
=======================

StereoMolGraph is designed to work seamlessly with `RDKit <https://www.rdkit.org/>`_, a popular cheminformatics toolkit.  
:ref:`Graph2RDMol`_ and :ref:`RDMol2Graph`_ functions allow for easy conversion between StereoMolGraph's graph representations and RDKit's :ref:`RDMol` objects.
only **Explicit Hydrogens**: StereoMolGraph uses only explicit hydrogens! 


Included Properties
--------------------
- Atomic connectivity (atoms and bonds).  
- Local stereochemistry (e.g., tetrahedral chirality).




Not included Properties
-------------------------
- Bond orders (single/double/triple)
- Formal atomic charges
- Unpaired electrons

Those are required they can be calculated using :ref:`stereomolgraph.algorithms.bond_orders.connectivity2bond`_