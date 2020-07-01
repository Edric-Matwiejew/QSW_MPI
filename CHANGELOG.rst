=========
Changelog
=========

**********************************
1.0.0 - 2020-06-26 Feature Release
**********************************

This is a major revision of QSW_MPI. The focus of this realse is the expansion of the simulation capabilities of QSW_MPI while focussing the scope of the package through removal of features which are better supported through pre-exisiting alernatives (specifically file I/O and visualisation).

Added
-----

* Generalised support for quantum stochsatic walks, including the non-moralising quantum stochastic walk through the :class:`~qsw_mpi.MPI.LQSW` and :class:`~qsw_mpi.MPI.GQSW` classes.
* Experimental support for sparse systems following the Gorini–Kossakowski–Sudarshan–Lindblad equation in its diagonalised form thorugh the :class:`~qsw_mpi.MPI.GKSL` class.
* Support for MPI-enabled parallel output to HDF5 using H5Py via the non-user accessible module :mod:`~qsw_mpi.parallel_io`.

* Addtional operator types including the cannonical Markov chain transition matrix, and those required for the demoralisation correction scheme.

Changed
-------

* All simulation types are now subclasses a generalised :class:`~qsw_mpi.MPI.walk` class. This breaks compatibility with version 0.0.1.
* :meth:`~qsw_mpi.MPI.walk.step` and :meth:`~qsw_mpi.MPI.walk.series` have been simplified, gathering of simulation results, or saving of the simulation results is now carried out through the :meth:`~qsw_mpi.MPI.walk.gather_result`, :meth:`~qsw_mpi.MPI.walk.gather_populations`, :meth:`~qsw_mpi.MPI.save_result` or :meth:`~qsw_mpi.MPI.save_populations`.

Removed
-------

* Removed visualisation module :mod:`~qsw_mpi.plot`. For basic visualisation direct use of Matplotlib and Networkx is recommended.
* Removed dedicated I/O module :mod:`~qsw_mpi.io`. For HDF5 file operations, direct use of H5Py is recommended.

**********************************
0.0.1 - 2020-03-05 Initial Release
**********************************

.. Note::
   This version supports only quantum stochastic walk simulation with locally interacting Lindblad operators.

