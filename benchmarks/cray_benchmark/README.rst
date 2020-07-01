==================
QSW Cray Benchmark
==================

These tests are intented to run in a distributed computing environment which uses the Slurm Workload Manager. Results obtained from these tests correspond to Figure 14 in "QSW_MPI: a framework for parallel simulation of quantum stochastic walks".

Additional Requirements
-----------------------

* Bash
* NetworkX
* Singularity

Setup Instructions
------------------

This benchmark makes use of Singularity container managment. To build the container, obtain an access token from https://cloud.sylabs.io and place it in this directory under the filename 'sylabs-token', then run 'remote_build.sh'. This builds the container (qsw_mpi.sif) on the Sylabs servers.

Then, after entering a valid account name in 'base.slurm', the benchmark is initialised via 'RUN_BENCHMARKS.sh'.

Pipeline
--------

#. 'cray_graph_gen.py' creates the test graph set:

  * Line graph, 5050 vertices.
  * Square lattice, 3844.
  * Erdos-Renyi graph with :math:`NLog(N)` edges, 2020 vertices.
  * Complete graph, 400 vertices.

#. 'slurm_gen.sh' generates slurm launch scripts suitable for measuring MPI only and MPI + OpenMP performance.
#. These jobs are submitted by the 'launch.sh' script present in the pure (MPI only) and hybrid (MPI + OpenMP) folder.
#. Simulation results and run-times are saved to pure/output and hybrid/output.
