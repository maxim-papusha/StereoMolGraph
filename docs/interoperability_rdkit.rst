RDKit Interoperability
=======================

StereoMolGraph is designed to work seamlessly with `RDKit <https://www.rdkit.org/>`_, a popular cheminformatics toolkit.  
:doc:`graph2rdmol </reference/graph2rdmol>` and :doc:`rdmol2graph </reference/rdmol2graph>` functions allow for easy conversion between StereoMolGraph's graph representations and RDKit's :ref:`RDMol` objects.
Only supports **explicit Hydrogens**!


Included Properties
--------------------
- Atomic connectivity (atoms and bonds).  
- Local stereochemistry (e.g., tetrahedral chirality).



Bond orders, Atomic charges and number of unpaired electrons can be calculated using :ref:`GenerateBondOrders <reference/stereomolgraph.algorithms.bond_orders.connectivity2bond>` algorithm.