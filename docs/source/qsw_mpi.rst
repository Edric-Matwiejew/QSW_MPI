Modules
=======

.. note::

   * :math:`N` is the dimensions of :math:`H, L` and :math:`\rho(t)`.
   * :math:`M` is the number of sources and sinks.
   * :math:`\tilde{N} = N + M`.

MPI
---

.. warning::
   The contents of this module will only operate correctly if called within an MPI instance.

.. automodule:: qsw_mpi.MPI
   :members:
   :undoc-members:
   :show-inheritance:

io
--

.. automodule:: qsw_mpi.io
   :members:
   :undoc-members:
   :show-inheritance:

measure
-------

.. automodule:: qsw_mpi.measure
   :members:
   :undoc-members:
   :show-inheritance:

operators
---------

.. automodule:: qsw_mpi.operators
   :members:
   :undoc-members:
   :show-inheritance:


.. _projection-anchor:

plot
----

:py:mod:`plot` provides basic visualization methods using `Matplotlib <https://matplotlib.org/>`_ as the plotting backend and `NetworkX <https://networkx.github.io>`_ to draw graphs. To use this module Matplotlib package should be imported:

.. code-block:: python

    import matplotlib.pyplot as plt

Methods :meth:`~qsw_mpi.plot.population_lines` and :meth:`~qsw_mpi.plot.coherence_lines` require a `Matplotlib axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ object as the first argument. This is instantiated as follows:

.. code-block:: python

    fig, ax = plt.subplot()


Methods :meth:`~qsw_mpi.plot.population_bars` and :meth:`~qsw_mpi.plot.coherence_bars` require a `Matplotlib axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ with a **3D projection**:

.. code-block:: python

    fig, ax = plt.subplot(projection = '3d')

Those unfamilar with the Matplotlib pyplot interface  are encouraged to consult the official `pyplot tutorial <https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html>`_.

.. automodule:: qsw_mpi.plot
   :members:
   :undoc-members:
   :show-inheritance:
