Installation
============
libNeST requires Python 3.9 or newer. We recommend installing it into a
virtual environment, which keeps it isolated from the system Python.

From PyPI
---------

.. code:: console

    $ python3 -m venv .venv
    $ source .venv/bin/activate
    $ python3 -m pip install libnest

From source
-----------
For development, clone the repository and install it in editable mode:

.. code:: console

    $ git clone https://github.com/danielpecak/libnest.git
    $ cd libnest
    $ python3 -m pip install -e ".[test]"     # library + pytest
    $ python3 -m pip install -e ".[docs]"     # library + Sphinx toolchain

Dependencies
------------
NumPy, SciPy, Matplotlib and pandas are installed automatically.

The optional package `py-WDATA <https://pypi.org/project/wdata/>`_ is needed
only by the example scripts that read WDATA simulation output
(``examples/wdata_slices.py``, ``examples/tools_example.py``):

.. code:: console

    $ python3 -m pip install wdata

Check the installation
----------------------

.. code:: console

    $ python3 -c "import libnest; print(libnest.__version__)"
    $ pytest                                  # from the repository root

If importing Matplotlib or pandas fails with
``numpy.core.multiarray failed to import``, the packages were built against a
different NumPy version. Install libNeST into a fresh virtual environment.

Next Steps
----------
See :ref:`tutorial` for examples.
