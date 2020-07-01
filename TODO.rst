====
Todo
====

In approximate order of importance:

   #. Provide examples of quantum walks utilising the :class:`~qsw_mpi.MPI.GKSL` class.
   #. Implement parallel-io operation using the Fortran HDF5 interface to remove dependancy on H5Py.
   #. Implement the :class:`~qsw_mpi.MPI.walk` class and subclasses :class:`~qsw_mpi.MPI.LQSW`, :class:`~qsw_mpi.MPI.GQSW` and :class:`~qsw_mpi.MPI.GKSL` in Fortran to make possible python indepedant walk simulation.
   #. Provide comprehensive documentation for the Fortran subroutines, compileable via Doxygen.
   #. Implement dense BLAS operations to increase the scope of effcient walk simulation to dense systems.
