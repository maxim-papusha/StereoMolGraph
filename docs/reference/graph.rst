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

.. autosummary::
   :toctree:  # Directory where stub files will be generated
   :recursive:       # Include submodules if needed
   
   stereomolgraph.MolGraph
   stereomolgraph.StereoMolGraph
   stereomolgraph.CondensedReactionGraph
   stereomolgraph.StereoCondensedReactionGraph