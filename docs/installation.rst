Installation
============

StereoMolGraph can be installed on Linux, MacOS, or Windows environment with **Python 3.10 or later**.

Install from PyPI
------------------

.. code-block:: bash

   pip install stereomolgraph

Installing from GitHub
----------------------

To install the latest development version directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/maxim-papusha/StereoMolGraph.git

Development Installation (Editable Mode)
----------------------------------------

If you plan to modify the source code, clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/maxim-papusha/StereoMolGraph.git
   cd StereoMolGraph
   pip install -e ".[dev]"
   

Now you can also run the tests with pytest:

.. code-block:: bash
    
    pytest
