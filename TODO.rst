====
Todo
====

#. Include the GKSL examples in the repository.
#. Create a conda environment file.
#. Make h5py optional for simulations that do not use HDF5.
#. Implement the walk classes in Fortran so simulations can run without Python.
#. Document the Fortran subroutines with Doxygen.
#. Add dense BLAS operations for simulations of dense systems.
#. Profile the full-interval GKSL cooling slowdown, measuring adjoint setup,
   norm estimation and Taylor products separately.
#. Compare blocking and overlapped MPI exchange across nodes for the Ising and
   shift jumps. Measure runtime and peak memory per rank to decide whether to
   reuse the communication graph or change the default.
#. Evaluate QuOp's HIP sparse matrix-vector and packing kernels on AMD GPUs.
   Compare device-aware MPI with host staging and validate the cooling
   trajectory against the CPU result.
