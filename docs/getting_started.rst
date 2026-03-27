Getting Started
===============

Install from PyPI
------------------

.. code-block:: bash

   pip install stereomolgraph


Creating a Molecular Graph
--------------------------

Build a :class:`~stereomolgraph.MolGraph` by specifying atoms and bonds:

.. code-block:: python

   from stereomolgraph import MolGraph

   # Create a water molecule (H-O-H) with explicit hydrogens
   g = MolGraph()
   g.add_atom(0, "O")
   g.add_atom(1, "H")
   g.add_atom(2, "H")
   g.add_bond(0, 1)
   g.add_bond(0, 2)

Adding Stereochemistry
----------------------

Use :class:`~stereomolgraph.StereoMolGraph` to include stereochemical information:

.. code-block:: python

   from stereomolgraph import StereoMolGraph
   from stereomolgraph.stereodescriptors import Tetrahedral

   g = StereoMolGraph()
   # central carbon
   g.add_atom(0, "C")
   g.add_atom(1, "F")
   g.add_atom(2, "Cl")
   g.add_atom(3, "Br")
   g.add_atom(4, "H")
   for i in range(1, 5):
       g.add_bond(0, i)

   # Assign tetrahedral stereo to atom 0
   g.add_stereo(0, Tetrahedral(neighbors=(1, 2, 3, 4), parity=1))

Comparing Molecules
-------------------

Graphs support equality and hashing out of the box:

.. code-block:: python

   g1 == g2          # isomorphism check
   hash(g1) == hash(g2)  # fast approximate comparison

See the :doc:`API Reference </reference/graph>` for the full list of graph operations, or
explore the :doc:`Menschutkin example notebook </examples/menschutkin_xyz>` for a more
complete workflow using 3D coordinates.
