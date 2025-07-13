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
   
   :maxdepth: 1
   :hidden:
   
   graphs/molgraph
   graphs/stereomolgraph
   graphs/condensedreactiongraph
   graphs/stereocondensedreactiongraph