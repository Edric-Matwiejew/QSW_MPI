Introduction
============

QSW_MPI is a python package developed for the modeling of quantum stochastic walks, a generalization of the classical random walk and continuous-time quantum walk in the Lindblad formalism. This model allows for the study of a wide range of Markovian open quantum systems subject to varying degrees of incoherent scattering.  Consisting of a python interface built on parallelized Fortran libraries utilizing sparse data structures; QSW_MPI is scalable to massively parallel computers, making possible the simulation of many thousands of graph vertices. QSW_MPI also provides explicit support for an extension of the standard quantum stochastic walk model to the study of non-Hermitian absorption and emission processes.

.. image:: images/animation.gif
    :alt: Quantum Stochastic Walk on a 2,3-balanced tree graph.
    :align: center
