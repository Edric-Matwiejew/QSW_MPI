======================
Desktop Benchmarks
======================

These tests are intended to gauge the performance and accuracy of QSW_MPI in a workstation-like (desktop) environment with respect to the pre-exisitng packages QSWalk.m and QSWalk.jl. The generated plots correspond to the following Figures in "QSW_MPI: a framework for parallel simulation of quantum stochastic walks":

* local_step: Figures 8, 11 and 12.
* global_step: Figures 9 and 13.
* series: Figures 10 and 15.


Addtional Requirements
----------------------

* Bash
* NetworkX
* Julia with QSWalk.jl
* Mathematica with Wolfram script support and QSWalk.m

Instructions
------------

These benchmarks are initialised via 'RUN_BENCHMARKS.sh'. This file contains the parameter defining the maximum run-time allowed for each test, which is 300s by default. The number of MPI processes used is defined by \*.sh scripts present in the local_step, global_step and series folders.

The results, timings and plots for each of the tests are stored in their respectuve subfolders. These are:

* local_step: Performance and accuracy of L-QSW simulation at a single time-point (QSWalk.m, QSWalk.jl and QSW_MPI).
* global_step: Performance and accuracy of G-QSW simulation at a single time-point (QSWalk.m, QSWalk.jl and QSW_MPI).
* series: Performance and accuracy of time-series L-QSW simulation (QSW_MPI only).
