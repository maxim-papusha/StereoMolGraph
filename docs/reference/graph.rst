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
   * - ``Molgraph``
     - ✓
     -
     -
     -
   * - ``CondensedReactionGraph``
     - ✓
     - ✓
     -
     -
   * - ``StereoMolGraph``
     - ✓
     -
     - ✓
     -
   * - ``StereoCondensedReactionGraph``
     - ✓
     - ✓
     - ✓
     - ✓


+-----------------------------+------------------+--------------------------+-----------------+-------------+
| Graph Class                 | Connectivity     | Connectivity             | Stereo          | Stereo      |
|                             |                  | Change                   |                 | Change      |
+=============================+==================+==========================+=================+=============+
| ``Molgraph``                | ✓                |                          |                 |             |
+-----------------------------+------------------+--------------------------+-----------------+-------------+
| ``CondensedReactionGraph``  | ✓                | ✓                        |                 |             |
+-----------------------------+------------------+--------------------------+-----------------+-------------+
| ``StereoMolGraph``          | ✓                |                          | ✓               |             |
+-----------------------------+------------------+--------------------------+-----------------+-------------+
| ``StereoCondensedReactionGraph`` | ✓             | ✓                       | ✓               | ✓           |
+-----------------------------+------------------+--------------------------+-----------------+-------------+


.. autoclass:: stereomolgraph.AtomId
   :no-index-entry:

.. autoclass:: stereomolgraph.Bond
   :no-index-entry:


.. toctree::
   :maxdepth: 1
   :caption: Graph Types:
   
   graphs/molgraph
   graphs/stereomolgraph
   graphs/condensedreactiongraph
   graphs/stereocondensedreactiongraph