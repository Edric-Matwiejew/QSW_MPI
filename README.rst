=======
QSW_MPI
=======

.. image:: https://readthedocs.org/projects/qsw-mpi/badge/?version=latest
    :target: https://qsw-mpi.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Overview
--------

QSW_MPI provides an parallel framework for quantum stochastic simulation. For an overview of package usage and the theoretical basis of quantum stochastic walks please consult the `documentation <https://qsw-mpi.readthedocs.io/en/latest/>`_.

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
    * Networkx (In order to run usage examples.)

Installation
------------

After cloning the repository enter 'QSW_MPI/src' and build the Fortran shared object libraries:

.. code-block:: bash

    make

After this the QSW_MPI package may be used by importing the 'QSW_MPI' folder to python's system path at runtime:

.. code-block:: python

    import sys
    sys.path.append('path_to/QSW_MPI')
    import qsw_mpi as qsw

Or, to install 'QSW_MPI' as normal, in the 'QSW_MPI' folder generate a distribution archive:

.. code-block:: python

    python3 setup.py sdist bdist_wheel

Enter the newly created 'QSW_MPI/dist' folder which should contain the archive 'qsw_mpi-0.0.1.tar.gz'. For with the QSW_MPI can be installed using pip3:

.. code-block:: bash

    pip3 install qsw_mpi-0.0.1.tar.gz

Usage
-----
A usage example is included in 'QSW_MPI/examples'. The is run by issuing the terminal command:

.. code-block:: bash

    mpiexec -N 2 python3 example.py


