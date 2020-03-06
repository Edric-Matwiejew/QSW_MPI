=======
QSW_MPI
=======

|docs| |doi|

.. |docs| image:: https://readthedocs.org/projects/qsw-mpi/badge/?version=latest
    :target: https://qsw-mpi.readthedocs.io/en/latest/?badge=latest

.. |doi| image:: https://zenodo.org/badge/205545419.svg
   :target: https://zenodo.org/badge/latestdoi/205545419

Overview
--------

QSW_MPI provides an parallel framework for quantum stochastic simulation. For an overview of package usage and the theoretical basis of quantum stochastic walks please consult the `documentation <https://qsw-mpi.readthedocs.io/en/latest/>`_, or preprint `article <https://arxiv.org/pdf/2003.02450.pdf>`_.

Requirements
------------
* GNU make
* GNU Fortran 5 or higher
* A functional MPI implementation
* HDF5
* Python 3.6.9 or higher with python packages packages:
    * mpi4py
    * NumPy
    * SciPy
    * H5py
    * Matplotlib
    * Networkx

Installation
------------

After cloning the repository enter 'QSW_MPI/src' and build the Fortran shared object libraries:

.. code-block::

    make

After this the QSW_MPI package may be used by importing the 'QSW_MPI' folder to python's system path at runtime:

.. code-block:: python

    import sys
    sys.path.append('path_to/QSW_MPI')
    import qsw_mpi as qsw

Or, to install 'QSW_MPI' as normal, in the 'QSW_MPI' folder generate a distribution archive:

.. code-block::

    python3 setup.py sdist bdist_wheel

Enter the newly created 'QSW_MPI/dist' folder which should contain the archive 'qsw_mpi-0.0.1.tar.gz'. For with the QSW_MPI can be installed using pip3:

.. code-block::

    pip3 install qsw_mpi-0.0.1.tar.gz

Documention
-----------

To obtain a local copy of the documentaion, with python package `Sphinx <http://www.sphinx-doc.org/en/master/>`_ and the Read the Docs Sphinx `Theme <https://sphinx-rtd-theme.readthedocs.io/en/stable/>`_ installed, enter `QSW\_MPI/docs` and build the documentation:

.. code-block::

    make html


Usage
-----
A usage example is included in 'QSW_MPI/examples'. The is run by issuing the terminal command:

.. code-block::

    mpiexec -N 2 python3 example.py

QSW_MPI Package Contents Overview
---------------------------------

Program Files
^^^^^^^^^^^^^
* qsw_mpi/
    * __init__.py - Python package initialization.
    * MPI.py - Parallel operations, quantum stochastic walk system creation and propagation.
    * operators.py - Creation of local quantum stochastic walk operators.
    * measure.py - Extract results from propagated walks.
    * io.py - Input and output of results.
    * plot.py - Basic visualization of results.

* src/
    * Makefile - Makefile to produce foperators and fMPI shared object libraries.
    * foperators.f90 - Source code for foperators shared object library.
    * fMPI.f90 - Source code for fMPI shared object library.
    * iso_precisions.f90 - Defines fortran precision types.
    * sparse.f90 - Sparse data representation and parallel BLAS operations.
    * one_norms.f90 - Parallel 1-norm estimation.
    * expm.f90 - Parallel calculation of the action of the matrix exponentional on a complex vector.
    * operatots.f90 - Creation of local and distributed quantum stochastic walk operators.

Other Files
^^^^^^^^^^^

* README.rst - QSW_MPI basic information.
* LICENSE - QSW_MPI license.
* setup.py - Configuration file used to generate a distribution archieve.
* MANIFEST.in - Additional files to include in the distribution archive.

* examples/
    * example.py - Usage example detailed in "QSW_MPI: A framework for parallel simulation of quantum stochastic walks".

* benchmark/
    * RUN_ME.sh - Benchmark automation bash script.
    * graphgen.py - Creation of test graph sets.
    * steps.py - Performs a quantum stochastic walks on the generated graphs.
    * plot_results.py - Plots time as a function of MPI processes.

* docs/
    * Makefile - Documentaion make script for Unix-like systems.
    * make.bat - Documenation build script for Windows systems.
    * requirements.txt - Requirements to build documentation of Read the Docs.
    * source/ - Documenation source files and images.
