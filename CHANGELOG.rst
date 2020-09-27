=========
Changelog
=========

****************************************************
1.0.1 - 2020-09-25 Bug Fix and Minor Feature Release
****************************************************

Added
-----

* All :class:`~qsw_mpi.MPI.walk` sub-classes will now initialise :math:`\rho(0)` as a maximally mixed state on passing `'mixed'` to :meth:`~qsw_mpi.MPI.walk.initial_state`.

Fixed
-----

* Fixed incorrect superoperator generation when using the :class:`~qsw_mpi.MPI.GQSW` class which occured when a global Lindblad operator had an empty row aligning with the first row of the local partition of the :math:`L^TL^* \otimes I_N` term in :math:`\tilde{\mathcal{L}}`.
* Edited the 'Theory' and 'Package Overview' sections of the package documentation to be in line with reviewer comments on the corresponding journal article.

Changed
-------

* Changed plot sizes and plot text sizes in response to reviewer feedback.
* Enabled plot legends in desktop benchmarks output.
* Desktop accuracy benchmarks now consider only :math:`\mathrm{max}(|\Delta\rho(t)|)`.

**********************************
1.0.0 - 2020-07-01 Feature Release
**********************************

This is a major revision of QSW_MPI. The focus of this release is the expansion of the simulation capabilities of QSW_MPI while focussing the scope of the package through the removal of features which are better supported through pre-existing alternatives (specifically file I/O and visualisation).

Added
-----

* Generalised support for quantum stochastic walks, including the non-moralising quantum stochastic walk through the :class:`~qsw_mpi.MPI.LQSW` and :class:`~qsw_mpi.MPI.GQSW` classes.
* Experimental support for sparse systems following the Gorini–Kossakowski–Sudarshan–Lindblad equation in its diagonalised form through the :class:`~qsw_mpi.MPI.GKSL` class.
* Support for MPI-enabled parallel output to HDF5 using H5Py via the non-user accessible module :mod:`~qsw_mpi.parallel_io`.

* Additional operator types including the canonical Markov chain transition matrix, and those required for the demoralisation correction scheme.

Changed
-------

* All simulation types are now subclasses a generalised :class:`~qsw_mpi.MPI.walk` class. This breaks compatibility with version 0.0.1.
* :meth:`~qsw_mpi.MPI.walk.step` and :meth:`~qsw_mpi.MPI.walk.series` have been simplified, gathering of simulation results, or saving of the simulation results is now carried out through the :meth:`~qsw_mpi.MPI.walk.gather_result`, :meth:`~qsw_mpi.MPI.walk.gather_populations`, :meth:`~qsw_mpi.MPI.save_result` or :meth:`~qsw_mpi.MPI.save_populations`.

Removed
-------

* Removed visualisation module :mod:`~qsw_mpi.plot`. For basic visualisation, direct use of Matplotlib and Networkx is recommended.
* Removed dedicated I/O module :mod:`~qsw_mpi.io`. For HDF5 file operations, direct use of H5Py is recommended.

**********************************
0.0.1 - 2020-03-05 Initial Release
**********************************

.. Note::
   This version supports only quantum stochastic walk simulation with locally interacting Lindblad operators.

