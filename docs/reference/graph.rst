Graphs
==================

.. list-table:: Graph Classes Comparison
   :widths: 30 20 20 15 15
   :header-rows: 1

   * - Graph Class
     - Connectivity
     - Connectivity Change
     - Stereo
     - Stereo Change
   * - :class:`~stereomolgraph.MolGraph`
     - ✓
     -
     -
     -
   * - :class:`~stereomolgraph.CondensedReactionGraph`
     - ✓
     - ✓
     -
     -
   * - :class:`~stereomolgraph.StereoMolGraph`
     - ✓
     -
     - ✓
     -
   * - :class:`~stereomolgraph.StereoCondensedReactionGraph`
     - ✓
     - ✓
     - ✓
     - ✓


.. toctree::
   :maxdepth: 2
   
   .. autoclass:: stereomolgraph.MolGraph
   :inherited-members:
   :show-inheritance:
   :member-order: bysource

  .. autoclass:: stereomolgraph.StereoMolGraph
   :inherited-members:
   :show-inheritance:
   :member-order: bysource

  .. autoclass:: stereomolgraph.CondensedReactionGraph
   :members:
   :inherited-members:
   :show-inheritance:
   :member-order: groupwise

  .. autoclass:: stereomolgraph.StereoCondensedReactionGraph
   :members:
   :inherited-members:
   :show-inheritance:
   :member-order: groupwise