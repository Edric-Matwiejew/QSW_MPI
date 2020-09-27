====
Todo
====

In approximate order of importance:

   #. Provide examples of quantum walks utilising the :class:`~qsw_mpi.MPI.GKSL` class.
   #. Create an environment file for Anaconda.
   #. Implement parallel-io operation using the Fortran HDF5 interface to remove the dependency on H5Py.
   #. Implement the :class:`~qsw_mpi.MPI.walk` class and subclasses :class:`~qsw_mpi.MPI.LQSW`, :class:`~qsw_mpi.MPI.GQSW` and :class:`~qsw_mpi.MPI.GKSL` in Fortran to make possible python indepedant walk simulation.
   #. Provide comprehensive documentation for the Fortran subroutines, compilable via Doxygen.
   #. Implement dense BLAS operations to increase the scope of efficient walk simulation to dense systems.
