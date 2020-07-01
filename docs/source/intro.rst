Introduction
============

QSW\_MPI is a python package developed for time-series simulation of continuous-time quantum stochastic walks. This model allows for the study of Markovian open quantum systems in the Lindblad formalism, including a generalisation of the continuous-time random walk and continuous-time quantum walk. Consisting of a python interface accessing parallelised Fortran libraries utilising sparse data structures, QSW\_MPI is scalable to massively parallel computers, which makes possible the simulation of a wide range of walk dynamics on directed and undirected graphs of arbitrary complexity. 

.. image:: graphics/animation.gif
    :alt: Quantum Stochastic Walk on a 2,3-balanced tree graph with a source and the centre vertex and sinks at the outer branches.
    :align: center
