#   QSW_MPI -  A package for parallel Quantum Stochastic Walk simulation.
#   Copyright (C) 2019 Edric Matwiejew
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import h5py
import qsw_mpi.operators as operators
import qsw_mpi.fMPI as fMPI
from qsw_mpi.operators import nm_pop_inv_map
import qsw_mpi.parallel_io as pio
from time import time
from warnings import warn

_MPI_COUNT_MAX = int(np.iinfo(np.int32).max)


def _nnz_offsets(comm, local_nnz):
    """
    Returns the zero-based offset of each rank's nonzero entries.
    The final offset is the total number of stored entries.
    """
    counts = np.empty(comm.size, dtype=np.int64)
    comm.Allgather(np.int64(local_nnz), counts)
    if np.any(counts < 0):
        raise MemoryError("A rank could not copy its assembled GKSL matrix")
    offsets = np.empty(comm.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    return offsets


def _gather_csr_array(comm, local, offsets, output, root):
    """
    Gathers a distributed CSR array to the root MPI process in blocks
    with 32-bit counts and displacements.
    """
    total = int(offsets[-1])
    for start in range(0, total, _MPI_COUNT_MAX):
        end = min(start + _MPI_COUNT_MAX, total)
        # MPI displacements are relative to the start of the current block.
        cuts = np.clip(offsets, start, end) - start
        counts = np.diff(cuts).astype(np.int32)
        local_start = max(0, start - int(offsets[comm.rank]))
        send = local[local_start:local_start + int(counts[comm.rank])]
        recv = (output[start:end], (counts, cuts[:-1].astype(np.int32))) if comm.rank == root else None
        comm.Gatherv(send, recv, root=root)


class walk(object):

    """
    The :class:`walk` class provides a framework for Lindbladian quantum walk simulation. This proceeds by creating a :class:`walk` sub-class containing a private method with the signature,

    .. code-block:: python

        def _construct_superoperator(self):

    which is responsible for defining:

    * self.SO_nnz (integer): The number of non-zeros in :math:`\\tilde{\mathcal{L}}`
    * self.SO_rows (integer): The global number of rows in :math:`\\tilde{\mathcal{L}}`
    * self.partition_table (integer array): A 1-based array describing the number of :math:`\\tilde{\mathcal{L}}` rows per MPI rank.

        :Example:
            | An evenly distributed :math:`\\tilde{\mathcal{L}}` with an MPI communicator of size of 4 and self.SO_nnz = 16 would have,
            |
            | self.partition_table = [1,5,9,13,17]
            |
            | where the number of rows at MPI rank i (self.SO_local_rows) is equal to self.partition_table[i + 1] - self.partition_table[i]

    * self.SO_local_rows (integer) the number of :math:`\\tilde{\mathcal{L}}` stored locally, :math:`N_\\text{local}`.
    * self.SO_n_row_starts (:math:`N_\\text{local} + 1` integer NumPy array): Local partition of the 1-based row index pointer array for :math:`\\tilde{\mathcal{L}}`, part of the CSR matrix format.
    * self.SO_col_indexes (:math:`N_\\text{local}` integer NumPy array): Local partition of the 1-based column index array for :math:`\\tilde{\mathcal{L}}`, part of the CSR matrix format.
    * self.SO_values (:math:`N_\\text{local}` complex NumPy array): Local partition of non-zero values corresponding to self.SO_col_indexes (in 1-based indexing) for :math:`\\tilde{\mathcal{L}}`, part of the CSR matrix format.

    * self.size (integer tuple): Dimensions of :math:`\\rho(t), H`, Lindblad and other non-vectorised operators.
    * self.MPI_communicator: MPI communicator object as defined by mpi4py.
    * self.rank (integer): MPI process rank


    :meth:`_construct_superoperator()` must be called in the '__init__' function of the :class:`walk` subclass, followed by a call to the `__init__` method of the parent class:

    .. code-block:: python

        super().__init__()

    By default, :meth:`step` and :meth:`series` use a blocking all-to-all
    exchange. Set ``overlap_communication`` to ``True`` on all MPI ranks to
    overlap the neighbour exchange with local matrix-vector multiplication.

    The following :class:`walk` sub-classes are defined in QSW_MPI:

    * :class:`LQSW`: QSWs with locally interacting Lindbad operators.
    * :class:`GQSW`: QSWs with globally interacting Lindblad operators (including non-moralising QSWs).
    * :class:`GKSL`: Quantum walks defined using the diagonal form of the GKSL equation.
    """

    def __init__(self):

        self.parallel_io = False
        self.overlap_communication = False
        self.local_result = None
        self._reconcile_mpi_communications()
        self._superoperator_one_norms()

    def _matricize_index(self, index, offset):
        """
        Takes a local index of the distributed :math:`\tilde{\\rho}(t)` and the offset of that index relative to its global position (self.partition_table[self.rank] - 1) and returns the corresponding indexes of :math:`\\rho(t)`.
        """
        j = np.int64(np.ceil((1 + index + offset)/self.size[0])) - 1
        i = np.int64(index + offset - j*self.size[0])

        return i, j

    def _check_dataset_name(self, File, dataset_name, action):
        """
        Checks if a dataset already exists in a HDF5 file, if so, '_' is append to the input name.
        """

        if (action == "w") or not (dataset_name in File):
            return dataset_name

        else:
            return dataset_name + '_'

    def _reconcile_mpi_communications(self):
        """
        Analyses the structure of :math:`\\tilde{\mathcal{L}}` to optimise MPI communication of :math:`\\tilde{\\rho}(t)` during sparse-matrix multiplication.
        """
        self.rec_time = time()
        self._kernel_row_starts = getattr(self, '_local_row_starts', self.SO_row_starts)

        def receive_send_indices(indices):
            try:
                self.SO_rhs_send_inds = indices.copy()
            except MemoryError:
                self.SO_rhs_send_inds = None

        (self.SO_local_col_inds, self.SO_num_rec_inds, self.SO_rec_disps,
         self.SO_num_send_inds, self.SO_send_disps) = fMPI.build_communications(
                self.SO_col_indexes, self.partition_table,
                self.MPI_communicator.py2f(), receive_send_indices)

        if self.MPI_communicator.allreduce(self.SO_rhs_send_inds is None, op=MPI.LOR):
            raise MemoryError("A rank could not copy its MPI send indices")
        self.rec_time = time() - self.rec_time

    def _superoperator_one_norms(self):
        """
        Calculates the 1-norm of the superoperator. :meth:`step` and
        :meth:`series` estimate the norms of higher powers when required.
        """
        self.one_norms_time = time()

        self.one_norms, p = fMPI.one_norm_series(
                self.SO_rows,
                self._kernel_row_starts,
                self.SO_col_indexes,
                self.SO_values,
                self.SO_num_rec_inds,
                self.SO_rec_disps,
                self.SO_num_send_inds,
                self.SO_send_disps,
                self.SO_local_col_inds,
                self.SO_rhs_send_inds,
                self.partition_table,
                self.MPI_communicator.py2f(), first_only=1)
        self.p = np.array(p, dtype=np.int32)

        self.one_norms_time = time() - self.one_norms_time

    def set_omega(self, omega):
        """
        Re-defines :math:`\omega` in the QSW master equation. Rebuilds the
        superoperator and its communication pattern and recalculates its 1-norm.

        :param omega: Interpolation parameter, :math:`\omega`
        :type omega: float
        """

        self.omega = omega
        self._construct_superoperator()
        self._reconcile_mpi_communications()
        self._superoperator_one_norms()

    def initial_state(self, state):
        """
        Sets the initial state :math:`\\rho(t_0)` and generates :math:`\\tilde{\\rho}(t_0)`.

        :param state: One of the following:

            * :math:`\\rho(t_0)` as a complex NumPy array.

            * :math:`p(t_0)`. An array describing the diagonal of :math:`\\rho(t_0)` which is then assumed to have no initial inter-vertex coherences.

            * 'mixed' - Equally distribute :math:`\\rho(t_0)` over all sites, including specified sources and sinks.
        """

        if not self.rho_set:

            if isinstance(state, str):

                if state == 'mixed':
                    self.rho = np.full((self.size[0]), np.float64(1)/np.float64(self.size[0]), dtype = np.complex128)
                    self.rho = np.diag(self.rho)
                    self.rho_set = True

                if not self.rho_set:
                    raise ValueError('Initial state name "' + state + '" not recognised by walk class.')

            else:

                if not isinstance(state, np.ndarray):
                    state = np.array(state, dtype = np.complex128)

                if len(state.shape) == 2:
                    self.rho = state
                    self.rho_set = True

                elif len(state.shape) == 1:
                    self.rho = np.diag(state)
                    self.rho_set = True

                if not self.rho_set:
                    raise ValueError('Unable to construct initial state from input array.')


        self.rho_v = fMPI.initial_state(
                self.SO_local_rows,
                self.rho,
                self.rank,
                self.partition_table,
                self.MPI_communicator.py2f())

    def initial_pure_state(self, amplitudes, root=0):
        r"""
        Set the distributed initial density operator to :math:`|\psi\rangle\langle\psi|`.

        Call this method on every rank of the walk's communicator with the same
        root. The root supplies a one-dimensional, unit-norm vector of complex
        amplitudes in the Hamiltonian basis. Other ranks may pass ``None``.
        For an LQSW with sources or sinks, a ket on the original graph is
        padded with zeros. A full ket may include amplitudes on the added sites.

        Each rank receives the ket and stores only its local slice of the
        column-stacked density operator in ``rho_v``.
        """
        comm = self.MPI_communicator
        if not isinstance(root, (int, np.integer)) or not 0 <= root < comm.size:
            raise ValueError("root must be a rank in the walk's communicator")

        n = self.size[0]
        error = None
        try:
            ket = np.empty(n, dtype=np.complex128)
            if self.rank == root:
                state = np.asarray(amplitudes)
                graph_size = n - getattr(self, 'pad', 0)
                if state.ndim != 1 or state.size not in (graph_size, n):
                    raise ValueError("Expected one amplitude per basis state "
                                     "(added source and sink sites may be omitted)")
                ket[:state.size] = state
                ket[state.size:] = 0
                if not np.isclose(np.linalg.norm(ket), 1.0, rtol=1e-6, atol=0):
                    raise ValueError("Pure-state amplitudes must have unit norm")
            local = np.empty(self.SO_local_rows, dtype=np.complex128)
        except Exception as exc:
            error = exc
        if comm.allreduce(error is not None, op=MPI.LOR):
            raise next(exc for exc in comm.allgather(error) if exc is not None)
        comm.Bcast(ket, root=root)

        column, row = divmod(int(self.partition_table[self.rank]) - 1, n)
        written = 0
        while written < local.size:
            count = min(n - row, local.size - written)
            np.multiply(ket[row:row + count], ket[column].conjugate(),
                        out=local[written:written + count])
            written += count
            column += 1
            row = 0

        self.rho_v = local
        self.rho = None
        self.rho_set = False
        self.local_result = None

    def step(self, t, precision = "dp"):
        """
        Calculate :math:`\\rho(t)`. This is done via an implementation of Algorithm 3.2 (without optional balancing or minimisation of the Frobenius norm) as detailed in:

         A. Al-Mohy and N. Higham, SIAM Journal on Scientific Computing 33, 488 (2011).

        :param t: time, :math:`t`.
        :type t: float

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision.
        :type precision: string, optional

        .. Note::
           To gather or save :math:`\\rho(t)` or :math:`p(t)`:

           * :meth:`~qsw_mpi.MPI.walk.gather_result`
           * :meth:`~qsw_mpi.MPI.walk.gather_populations`
           * :meth:`~qsw_mpi.MPI.walk.save_result`
           * :meth:`~qsw_mpi.MPI.walk.save_populations`
        """
        self.step_time = time()

        self.local_result = fMPI.step(
                self.SO_rows,
                self.SO_local_rows,
                self._kernel_row_starts,
                self.SO_col_indexes,
                self.SO_values,
                self.SO_num_rec_inds,
                self.SO_rec_disps,
                self.SO_num_send_inds,
                self.SO_send_disps,
                self.SO_local_col_inds,
                self.SO_rhs_send_inds,
                t,
                self.rho_v,
                self.partition_table,
                self.p,
                self.one_norms,
                self.MPI_communicator.py2f(),
                precision, overlap=self.overlap_communication)

        self.step_time = time() - self.step_time

    def series(self, t1, tq, steps, precision = 'dp'):
        """
        Calculate :math:`\\rho(t)` over :math:`[t_1, t_q]` (where :math:`t_q > t_1`) at :math:`q` evenly spaced steps with :math:`\Delta t = \\frac{t_q - t_1}{q}`. This is done via an implementation of Algorithm 5.2 (without optional balancing or minimisation of the Frobenius norm) as detailed in:

         A. Al-Mohy and N. Higham, SIAM Journal on Scientific Computing 33, 488 (2011).

        :param t1: Starting time, :math:`t_1`.
        :type t1: float

        :param tq: Ending time, :math:`t_q`.
        :type tq: float

        :param steps: :math:`q`.
        :type steps: integer

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision. Defaults to "dp".
        :type precision: string, optional

        .. Note::
           To gather or save :math:`(\\rho(t_0),...,\\rho(t_q))` or :math:`(p(t_0),...,p(t_0))`:

           * :meth:`~qsw_mpi.MPI.walk.gather_result`
           * :meth:`~qsw_mpi.MPI.walk.gather_populations`
           * :meth:`~qsw_mpi.MPI.walk.save_result`
           * :meth:`~qsw_mpi.MPI.walk.save_populations`
        """
        self.series_time = time()

        self.local_result = fMPI.series(
                self.SO_rows,
                self.SO_local_rows,
                self._kernel_row_starts,
                self.SO_col_indexes,
                self.SO_values,
                self.SO_num_rec_inds,
                self.SO_rec_disps,
                self.SO_num_send_inds,
                self.SO_send_disps,
                self.SO_local_col_inds,
                self.SO_rhs_send_inds,
                t1,
                tq,
                self.rho_v,
                steps,
                self.partition_table,
                self.p,
                self.one_norms,
                self.MPI_communicator.py2f(),
                precision, overlap=self.overlap_communication)

        self.series_time = time() - self.series_time

    def gather_result(self, root = 0):
        """
        Gathers :math:`\\rho(t)` or :math:`(\\rho(t_1),...,\\rho(t_q))` following execution of :meth:`qsw_mpi.MPI.walk.step` or :meth:`qsw_mpi.MPI.walk.series` at a specified MPI rank. `None` is returned elsewhere.

        :param root: MPI process rank at which to gather :math:`\\rho(t)` or :math:`(\\rho(t_1),...,\\rho(t_q))`.
        :type root: integer

        :returns: :math:`\\tilde{N} \\times \\tilde{N}` complex NumPy array at root, otherwise `None`.
        """

        if self.local_result is not None:

            if self.rank == root:
                output_size = self.size
            else:
                output_size = 1

            if self.local_result.ndim == 1:

                rhot = fMPI.gather_step(
                        self.local_result,
                        self.partition_table,
                        root,
                        self.MPI_communicator.py2f(),
                        output_size).T

            else:

                rhot = fMPI.gather_series(
                            self.local_result,
                            self.partition_table,
                            root,
                            self.MPI_communicator.py2f(),
                            output_size).T

            return rhot

    def _get_local_populations(self):
        """
        Returns portion of :\\vec{p}: corresponding to the local partition of :math:`\\tilde{\mathcal{L}}`.
        """

        stride = self.size[0] + 1
        offset = int(self.partition_table[self.rank]) - 1
        return np.ascontiguousarray(self.local_result[(-offset) % stride::stride])

    def gather_populations(self, root = 0):
        """
        Gathers :math:`p(t)` or :math:`((p(t_0),...,p(t_q))` following execution of :meth:`qsw_mpi.MPI.walk.step` or :meth:`qsw_mpi.MPI.walk.series` at a specified MPI rank. `None`` is returned elsewhere.

        :param root: MPI process rank at which to gather :math:`p(t)` or :math:`(p(t_1),...,p(t_q))`.
        :type root: integer

        :returns: :math:`p(t)` or :math:`(p(t_0),...,p(t_q))`.
        :rtype: :math:`N`, or :math:`(q + 1) \\times N`, complex NumPy array.
        """
        if self.local_result.ndim == 1:

            local_populations = self._get_local_populations()

            send_counts = np.array(self.MPI_communicator.gather(len(local_populations),root))

            if self.rank == root:
                populations = np.empty(np.sum(send_counts), dtype = np.complex128)
            else:
                populations = None

            self.MPI_communicator.Gatherv(local_populations, (populations, send_counts), root)

            return populations

        else:

            local_populations = self._get_local_populations()

            send_counts = np.array(self.MPI_communicator.gather(local_populations.shape[0] * local_populations.shape[1], root))

            disps = np.empty(self.flock, dtype = int)
            disps[:] = 0

            if self.rank == root:
                for i in range(self.flock - 1):
                    disps[i + 1] = disps[i] + send_counts[i]

            if self.rank == root:
                populations = np.empty((self.size[0], self.local_result.shape[1]), dtype = np.complex128)
            else:
                populations = None

            self.MPI_communicator.Gatherv(local_populations, [populations, send_counts, disps, MPI.DOUBLE_COMPLEX], root)

            return populations

    def gather_superoperator(self, root = 0):
        """
        Gathers :math:`\\tilde{\mathcal{L}}` to the root MPI process as a SciPy CSR matrix.

        :param root: MPI process rank at which to gather :math:`\\tilde{\mathcal{L}}`.
        :type root: integer

        :returns: :math:`\\tilde{\mathcal{L}}`
        :rtype: :math:`\\tilde{N}^2 \\times \\tilde{N}^2` SciPY CSR matrix
        """

        comm = self.MPI_communicator
        offsets = _nnz_offsets(comm, len(self.SO_values))
        row_offsets = self.partition_table.astype(np.int64) - 1
        # Include the final row pointer from the last rank.
        row_offsets[-1] += 1
        row_count = row_offsets[self.rank + 1] - row_offsets[self.rank]
        values = columns = pointers = None
        allocation_ok = True
        if self.rank == root:
            try:
                values = np.empty(offsets[-1], dtype=np.complex128)
                columns = np.empty(offsets[-1], dtype=np.int32)
                pointers = np.empty(row_offsets[-1], dtype=np.int64)
            except MemoryError:
                allocation_ok = False
        if not comm.bcast(allocation_ok, root=root):
            raise MemoryError("The gather root could not allocate the global CSR matrix")

        _gather_csr_array(comm, self.SO_values, offsets, values, root)
        _gather_csr_array(comm, self.SO_col_indexes, offsets, columns, root)
        _gather_csr_array(comm, self.SO_row_starts[:row_count].astype(np.int64, copy=False),
                          row_offsets, pointers, root)
        if self.rank == root:
            columns -= 1
            pointers -= 1
            return sp.csr_matrix((values, columns, pointers), shape=(self.SO_rows, self.SO_rows))

    def save_result(self, filename, walkname, action = 'a'):
        """
        Saves :math:`\\rho(t)` or :math:`\\vec{\\rho}(t)` to a HDF5 file.

        :param filename: File basename.
        :type filename: string

        :param walkname: Group/dataset name for :math:`\\rho(t)` or :math:`\\vec{\\rho}(t)`
        :type walkname: string

        :param action: File access type, append ('a'), overwrite ('w').
        :param action: string
        """
        if self.parallel_io:
            pio.save_result(self, filename, walkname, action)
        else:
            rhot = self.gather_result()

            if self.rank == 0:
                with h5py.File(filename + '.h5', action) as f:
                    f.create_dataset(self._check_dataset_name(f, walkname, action),
                                     data = rhot, dtype = pio.complex_dtype)

    def save_populations(self, filename, popname, action = 'a'):
        """
        Saves :math:`p(t)` or :math:`(p(t_1),...,p(t_q))` to a HDF5 file.

        :param filename: File basename.
        :type filename: string

        :param walkname: Group/dataset name for :math:`p(t)` or :math:`(p(t),...,p(t_q))`.
        :type walkname: string

        :param action: File access type, append ('a'), overwrite ('w').
        :param action: string
        """
        if self.parallel_io:
            pio.save_populations(self, filename, popname, action)
        else:

            populations = self.gather_populations()

            if self.rank == 0:
                with h5py.File(filename + '.h5', action) as f:
                    f.create_dataset(popname, data = populations, dtype = pio.complex_dtype)

    def save_superoperator(self, filename, operatorname, action = 'a'):
        """
        Save :math:`\\tilde{\mathcal{L}}` to a HDF5 file as a CSR matrix.

        :param filename: File basename
        :type filename: string

        :param operatorname: Group/dataset name for :math:`\\tilde{\mathcal{L}}`.
        :type operatorname: string

        :param action: File access type, append ('a'), overwrite ('w').
        :type action: string
        """
        if self.parallel_io:
            pio.save_superoperator(self, filename, operatorname, action)
        else:

            M = self.gather_superoperator()

            if self.rank == 0:
                with h5py.File(filename + '.h5', action) as f:
                    f.create_dataset(self._check_dataset_name(f, operatorname + '/indptr', action), data = M.indptr)
                    f.create_dataset(self._check_dataset_name(f, operatorname + '/indices', action), data = M.indices)
                    f.create_dataset(self._check_dataset_name(f, operatorname + '/data', action),
                                     data = M.data, dtype = pio.complex_dtype)
                    rows = f[operatorname + '/indptr'].shape[0] - 1
                    f[operatorname].attrs['dimensions'] = (rows, rows)

class LQSW(walk):
    """
    Defines an L-QSW, a QSW with locally interacting Lindblad operators, with or without appended sources and sinks (see Equation :eq:`eq:qsw` or :eq:`eq:qsw_ss`).

    :param omega: Interpolation parameter, :math:`\omega`
    :type omega: float

    :param H: Hamiltonian, :math:`H`
    :type H: :math:`N \\times N` SciPy CSR

    :param L: :math:`M_L = \sum_k L_k` where the Lindblad operators are defined as in Equation :eq:`eq:lindblad`.
    :type L: :math:`N \\times N` SciPy CSR

    :param MPI_communicator: MPI communicator object provided by mpi4py.

    :param sources: Contains a list of absorption sites and the corresponding absorption rates.
    :type sources: tuple ([integers],[floats]), optional

    :param sinks: Contains a list of the emission sites and the corresponding absorption rates.
    :type sinks: tuple ([integers],[floats]), optional

    :return: A :class:`walk` object with a :math:`\\tilde{\mathcal{L}}` distributed over an MPI communicator.
    """
    def __init__(
            self,
            omega,
            H,
            L,
            MPI_communicator,
            sources = (None, None),
            sinks = (None, None)):

        self.MPI_communicator = MPI_communicator
        self.flock = self.MPI_communicator.Get_size()
        self.rank = self.MPI_communicator.Get_rank()
        self.omega = omega

        self.H = H
        operators.check_indices(self.H)

        self.L = L
        operators.check_indices(self.L)

        self.source_sites = sources[0]
        self.source_rates = sources[1]

        self.sink_sites = sinks[0]
        self.sink_rates = sinks[1]

        if (self.source_sites is not None) and (self.sink_sites is not None):

            self.pad = len(self.source_sites) + len(self.sink_sites)
            self.n_sources = len(self.source_sites)
            self.n_sinks = len(self.sink_sites)

            source_idx = np.argsort(sources[0])
            self.source_sites = np.array(sources[0])[source_idx]
            self.source_rates = np.array(sources[1])[source_idx]

            sink_idx = np.argsort(sinks[0])
            self.sink_sites = np.array(sinks[0])[sink_idx]
            self.sink_rates = np.array(sinks[1])[sink_idx]

        elif self.source_sites is not None:

            self.pad = len(self.source_sites)
            self.n_sources = len(self.source_sites)
            self.n_sinks = 0
            self.sink_sites = [-1]
            self.sink_rates = [0]

            source_idx = np.argsort(sources[0])
            self.source_sites = np.array(sources[0])[source_idx]
            self.source_rates = np.array(sources[1])[source_idx]

            sink_idx = np.argsort(sinks[0])
            self.sink_sites = np.array(sinks[0])[sink_idx]
            self.sink_rates = np.array(sinks[1])[sink_idx]

        elif self.sink_sites is not None:

            self.pad = len(self.sink_sites)
            self.n_sinks = len(self.sink_sites)
            self.n_sources = 0
            self.source_sites = [-1]
            self.source_rates = [0]

        else:

            self.pad = 0
            self.n_sources = 0
            self.n_sinks = 0
            self.source_sites = [-1]
            self.source_rates = [0]
            self.sink_sites = [-1]
            self.sink_rates = [0]

        self.size = (H.shape[0] + self.pad, H.shape[1] + self.pad)

        if self.source_sites[-1] > H.shape[0] - 1:
            raise IndexError("Maximum source site index not in range of H.")

        if self.sink_sites[-1] > H.shape[0] - 1:
            raise IndexError("Maximum sink site index not in range of H.")

        self._construct_superoperator()

        super().__init__()

    def _construct_superoperator(self):
        """
        Constructs a distributed vectorised L-QSW superoperator, :math:`\\tilde{\mathcal{L}}`.
        """

        if (self.omega < 0) or (self.omega > 1):
            raise ValueError("Parameter omega must satify 0 <= omega =< 1.")

        self.construct_time = time()

        self.SO_nnz, self.SO_rows, self.partition_table = fMPI.super_operator_extent(
                self.omega,
                self.H.indptr,
                self.H.indices,
                self.H.data,
                self.L.indptr,
                self.L.indices,
                self.L.data,
                self.source_sites,
                self.source_rates,
                self.sink_sites,
                self.sink_rates,
                self.flock,
                self.MPI_communicator.py2f())

        self.SO_local_rows = self.partition_table[self.rank + 1] - self.partition_table[self.rank]
        self.SO_n_row_starts = np.max([2, self.SO_local_rows + 1])

        self.SO_row_starts, self.SO_col_indexes, self.SO_values = fMPI.super_operator(
                self.omega,
                self.H.indptr,
                self.H.indices,
                self.H.data,
                self.L.indptr,
                self.L.indices,
                self.L.data,
                self.source_sites,
                self.source_rates,
                self.sink_sites,
                self.sink_rates,
                self.SO_nnz,
                self.SO_n_row_starts,
                self.MPI_communicator.py2f())

        self.construct_time = time() - self.construct_time

    def initial_state(self, state):
        """
        Sets the initial state, :math:`\\rho(t_0)`. For an :class:`~qsw_mpi.MPI.LQSW`, this method accepts an additional keyword argument.

        :param state: One of the following:

            * :math:`\\rho(t_0)` as a complex NumPy array.

            * :math:`p(t_0)`. An array describing the diagonal of :math:`\\rho(t_0)` which is then assumed to have no initial inter-vertex coherences.

            * 'sources' - Equally distribute :math:`\\rho(t_0)` over all specified sources.

            * 'mixed' - Equally distribute :math:`\\rho(t_0)` over all sites, including specified sources and sinks.
        """

        self.rho_set = False

        if isinstance(state, str):

            if state == 'sources':
                self.rho = np.zeros((self.size[0]), dtype = np.complex128)
                self.rho[-(self.n_sinks + self.n_sources):-self.n_sinks] = np.float64(1)/np.float64(self.n_sources)
                self.rho = np.diag(self.rho)
                self.rho_set = True

        elif not isinstance(state, str):

            if not isinstance(state, np.ndarray):
                state = np.array(state, dtype = np.complex128)

            if len(state.shape) == 2:
                self.rho = np.concatenate((state, np.zeros((state.shape[0], self.pad))), axis = 1)
                self.rho = np.concatenate((self.rho, np.zeros((self.pad, self.rho.shape[1]))), axis = 0)
                self.rho_set = True

            elif len(state.shape) == 1:
                self.rho = np.diag(np.concatenate((state, np.zeros((self.pad,))), axis = 0))
                self.rho_set = True

            if not self.rho_set:
                raise ValueError('Unable to construct initial state from input array.')

        super().initial_state(state)

class GKSL(walk):
    """
    :class:`~qsw_mpi.MPI.GKSL` provides for the simulation of systems obeying the GKSL master equation in its diagonalised form (see Equation :eq:`eq:gksl`).

    ``SO_row_starts`` contains one-based, 64-bit global row pointers.
    Local row pointers, column indices and MPI counts are 32-bit. To allow
    for the final row pointer, the global row count and the number of stored
    entries per rank must each be less than :math:`2^{31}-1`. These limits
    also apply to the adjoint used in norm estimation.

    :param H: Hamiltonian, :math:`H`.
    :type H: :math:`N \\times N` complex SciPy CSR

    :param Ls: Generalised Lindblad operators :math:`L_k`.
    :type Ls: :math:`N \\times N` complex SciPy CSR matrix

    :param taus: Coefficients :math:`\\tau_k` of the Lindblad operators :math:`L_k` (see Equation :eq:`eq:KL_eq`). If `None`, :math:`\\tau_k = 1`.
    :type taus: float array, optional

    :param MPI_communicator: MPI communicator object provided by mpi4py.

    :return: A :class:`walk` object with :math:`\\tilde{\mathcal{L}}` distributed over an MPI communicator.
    """
    def __init__(
            self,
            H,
            Ls,
            MPI_communicator,
            taus = None):

        self.MPI_communicator = MPI_communicator
        self.flock = self.MPI_communicator.Get_size()
        self.rank = self.MPI_communicator.Get_rank()

        self.size = H.shape

        if taus is None:
            self.taus = np.ones(len(Ls))
        else:
            self.taus = taus

        self.H = H
        operators.check_indices(self.H)

        self.Ls = Ls

        for L in Ls:
            operators.check_indices(L)

        self._construct_superoperator()

        super().__init__()

    def _construct_superoperator(self):
        """
        Construct distributed :math:`\\tilde{\mathcal{L}}`.
        """
        self.construct_time = time()

        if self.size[0]**2 >= _MPI_COUNT_MAX:
            raise OverflowError("GKSL row count and its one-based sentinel must fit in int32")
        if (sum(L.nnz for L in self.Ls) >= _MPI_COUNT_MAX
                or len(self.Ls) * (self.size[0] + 1) > _MPI_COUNT_MAX):
            raise OverflowError("Packed Lindblad inputs exceed the int32 interface")

        # Pack Ls
        self.L_data = self.Ls[0].data
        self.L_indices = self.Ls[0].indices
        self.L_indptr = self.Ls[0].indptr
        self.L_nnz = [self.Ls[0].nnz]

        if len(self.Ls) > 1:
            for i in range(1,len(self.Ls)):
                self.L_data = np.concatenate((self.L_data, self.Ls[i].data))
                self.L_indices = np.concatenate((self.L_indices, self.Ls[i].indices))
                self.L_indptr = np.concatenate((self.L_indptr, self.Ls[i].indptr))
                self.L_nnz.append(self.Ls[i].nnz)

        def receive_csr(ptr, cols, vals):
            # Copy the CSR arrays before the Fortran routine deallocates them.
            try:
                self._local_row_starts = ptr.copy() if ptr.size > 1 else np.repeat(ptr, 2)
                self.SO_row_starts = self._local_row_starts.astype(np.int64)
                self.SO_col_indexes = cols.copy()
                self.SO_values = vals.copy()
            except MemoryError:
                self._local_row_starts = None

        self.SO_nnz, self.SO_rows, self.partition_table = fMPI.build_gksl(
                self.taus,
                self.H.indptr,
                self.H.indices,
                self.H.data,
                self.L_nnz,
                self.L_indptr,
                self.L_indices,
                self.L_data,
                self.flock,
                self.MPI_communicator.py2f(),
                receive_csr)

        copied_nnz = self.SO_nnz if self._local_row_starts is not None else -1
        offsets = _nnz_offsets(self.MPI_communicator, copied_nnz)
        self.SO_row_starts += offsets[self.rank]
        self.SO_local_rows = self.partition_table[self.rank + 1] - self.partition_table[self.rank]
        self.SO_n_row_starts = self._local_row_starts.size

        self.construct_time = time() - self.construct_time

    def set_omega(self, omega):
        warnings.warn("set_omega() not supported by GKSL class, superoperator unchanged.")

    def initial_state(self, state):

        """
        Sets the initial state :math:`\\rho(t_0)` and generates :math:`\\tilde{\\rho}(t_0)`.

        :param state: One of the following:

            * :math:`\\rho(t_0)` as a complex NumPy array.

            * :math:`p(t_0)`. An array describing the diagonal of :math:`\\rho(t_0)` which is then assumed to have no initial inter-vertex coherences.

            * 'mixed' - Equally distribute :math:`\\rho(t_0)` over all sites, including specified sources and sinks.
        """

        self.rho_set = False

        super().initial_state(state)

class GQSW(walk):
    """
    :class:`~qsw_mpi.MPI.GQSW` provides for the simulation of G-QSWs and NM-G-QSWs (see Equations :eq:`eq:L_global` or :eq:`eq:nm_gqsw`).

    :param omega: Interpolation parameter, :math:`\omega`.
    :type omega: float

    :param H: Hamiltonian, :math:`H`.
    :type H: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR

    :param Ls: Generalised Lindblad operators :math:`L_k`.
    :type Ls: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR matrix

    :param MPI_communicator: MPI communicator object provided by mpi4py.

    :param H_loc: NM-G-QSW only. Hamiltonian acting in the subspaces of :math:`V^D`, produced by :meth:`~qsw_mpi.operators.nm_H_loc` (see Equation :eq:`eq:nm_gqsw`).
    :type H_loc: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR, optional

    :param vsets: NM-G-QSW only. Set of vertex subspaces :math:`V^D` created by :meth:`~qsw_mpi.operators.nm_vsets`, see (Equation :eq:`eq:nm_gqsw`).
    :type vsets: integer lists, optional

    :return: A :class:`walk` object with a :math:`\\tilde{\mathcal{L}}` distributed over an MPI communicator.
    """
    def __init__(
            self,
            omega,
            H,
            Ls,
            MPI_communicator,
            H_loc = None,
            vsets = None):

        self.MPI_communicator = MPI_communicator
        self.flock = self.MPI_communicator.Get_size()
        self.rank = self.MPI_communicator.Get_rank()
        self.omega = omega

        self.size = H.shape

        self.H = H
        operators.check_indices(self.H)

        self.Ls = Ls

        for L in Ls:
            operators.check_indices(L)

        if H_loc is not None:
            self.H_loc = H_loc
            operators.check_indices(self.H_loc)
            self.H_loc_indptr = H_loc.indptr
            self.H_loc_indices = H_loc.indices
            self.H_loc_data = H_loc.data
        else:
            self.H_loc_indptr = np.zeros(H.shape[0] + 1, dtype = np.int32)
            self.H_loc_indices = np.zeros(1, dtype = np.int32)
            self.H_loc_data = np.zeros(1, dtype = np.complex128)

        self.vsets = vsets

        self._construct_superoperator()

        super().__init__()

    def _construct_superoperator(self):
        """
        Construct distributed :math:`\\tilde{\mathcal{L}}` for a G-QSW or NM-G_QW.
        """

        if (self.omega < 0) or (self.omega > 1):
            raise ValueError("Parameter omega must satify 0 <= omega =< 1.")

        self.construct_time = time()

        # Pack Ls
        self.L_data = self.Ls[0].data
        self.L_indices = self.Ls[0].indices
        self.L_indptr = self.Ls[0].indptr
        self.L_nnz = [self.Ls[0].nnz]

        if len(self.Ls) > 1:
            for i in range(1,len(self.Ls)):
                self.L_data = np.concatenate((self.L_data, self.Ls[i].data))
                self.L_indices = np.concatenate((self.L_indices, self.Ls[i].indices))
                self.L_indptr = np.concatenate((self.L_indptr, self.Ls[i].indptr))
                self.L_nnz.append(self.Ls[i].nnz)

        self.SO_nnz, self.SO_rows, self.partition_table = fMPI.demoralized_super_operator_extent(
                self.omega,
                self.H.indptr,
                self.H.indices,
                self.H.data,
                self.H_loc_indptr,
                self.H_loc_indices,
                self.H_loc_data,
                self.L_nnz,
                self.L_indptr,
                self.L_indices,
                self.L_data,
                self.flock,
                self.MPI_communicator.py2f())

        self.SO_local_rows = self.partition_table[self.rank + 1] - self.partition_table[self.rank]
        self.SO_n_row_starts = np.max([2, self.SO_local_rows + 1])

        self.SO_row_starts, self.SO_col_indexes, self.SO_values = fMPI.demoralized_super_operator(
                self.omega,
                self.H.indptr,
                self.H.indices,
                self.H.data,
                self.H_loc_indptr,
                self.H_loc_indices,
                self.H_loc_data,
                self.L_nnz,
                self.L_indptr,
                self.L_indices,
                self.L_data,
                self.partition_table,
                self.MPI_communicator.py2f(),
                self.SO_nnz,
                self.SO_n_row_starts)

        self.construct_time = time() - self.construct_time

    def initial_state(self, state):

        """
        Sets the initial state :math:`\\rho(t_0)` and generates :math:`\\tilde{\\rho}(t_0)`.

        :param state: One of the following:

            * :math:`\\rho(t_0)` as a complex NumPy array.

            * :math:`p(t_0)`. An array describing the diagonal of :math:`\\rho(t_0)` which is then assumed to have no initial inter-vertex coherences.

            * 'mixed' - Equally distribute :math:`\\rho(t_0)` over all sites, including specified sources and sinks.
        """

        self.rho_set = False

        super().initial_state(state)

    def gather_populations(self, root = 0):
        """
        Return :math:`p(t)` or :math:`(p(t_0),...,p(t_q))` following a call to :meth:`~qsw_mpi.MPI.walk.step` or :meth:`~qsw_mpi.MPI.walk.series`. In the case of a G-QSW, :math:`p(t)` corresponds to the originating graph :math:`G` (see Equation :eq:`eq:nm_gqsw`).

        :param root: MPI process rank at which to gather :math:`p(t)` or :math:`(p(t_0),...,p(t_q))`.
        :type root: integer, optional.

        :returns: :math:`p(t)` or :math:`(p(t_0),...,p(t_q))`.
        :rtype: :math:`N`, or :math:`(q + 1) \\times N`, complex NumPy array.
        """

        populations = super().gather_populations(root)

        if self.vsets == None:
            return populations

        if self.rank == root:
            return nm_pop_inv_map(populations, self.vsets)

    def save_vsets(self, filename, vsetsname, action = 'a'):
        """
        Save :math:`V^D` (see Equation :eq:`eq:nm_gqsw`) created by :meth:`~qsw_mpi.operators.nm_vsets` to a HDF5 file.

        :param filename: File basename.
        :type filename: string

        :param vsetname: Group/dataset name under which to save :math:`V^D`.
        :type vsetname: string

        :param action: File access type, append ('a'), overwrite ('w').
        :type action: string
        """

        if self.rank == 0:

            f = h5py.File(filename + '.h5', action)

            savename = self._check_dataset_name(f, vsetsname, action)

            for i, v in enumerate(self.vsets):

                dset = f.create_dataset(
                        vsetsname + '/' + str(i),
                        (len(v),),
                        dtype = np.int)

                dset[:] = np.array(v, dtype = np.int)

            f.close()

