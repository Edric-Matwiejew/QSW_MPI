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
import qsw_mpi.io as io
import time
from warnings import warn

def load_walk(
        omega,
        filename,
        MPI_communicator):
    """ Initialize a previously created :class:`walk` using a .qsw file.

    :param omega: :math:`\omega`
    :type omega: float

    :param filename: Target .qsw file to open.
    :type filename: string

    :param MPI_communicator: MPI communicator

    :return: A :class:`walk` object initialized using the input systems's :math:`H`, :math:`L`, sources, and sinks data. The :class:`walk.file` parameter is set to the input file.
    """
    H, L, sources, sinks = load_system(filename, MPI_communicator)

    walk_out = walk(omega, H, L, MPI_communicator, sources, sinks)
    walk_out.File(filename, action = "a")

    return walk_out

class walk(object):
    """The :class:`walk` object handels the creation of distributed :math:`\mathcal{L}`, :math:`\\rho(0)` and system propagation in an MPI environment. It also enables the saving of results and input data to a .qsw file. This object and its containing methods will only function correctly if called within an MPI instance, as shown in the examples.

    :param omega: :math:`\omega`
    :type omega: float

    :param H: :math:`H`
    :type H: (N, N), complex

    :param L: :math:`L`
    :type L: (N, N), complex

    :param MPI_communicator: MPI communicator

    :param sources: Contains a list of  absorption sites and the corresponding absorption rates.
    :type sources: ([M], [M]), (integer, float), tuple, optional

    :param sinks: Contains a list of the emission sites and the corresponding absorption rates.
    :type sinks: ([M], [M]), (integer, float), tuple, optional

    :return: A :class:`walk` object with a :math:`\mathcal{L}` distributed over an MPI communicator.
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

        self.filename = None

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

        self.construct_superoperator()

    def construct_superoperator(self):
        """A distributed :math:`\mathcal{L}` is constructed on :class:`walk` initialization\
                or on calling :py:meth:`~freeqsw.MPI.walk.set_omega`. One-norms \
                of the :math:`\mathcal{L}` power series are calculated for using in \
                deriving the optimal matrix exponentiation parameters when calling \
                :py:meth:`~freeqsw.MPI.walk.step` and :py:meth:`~freeqsw.MPI.walk.series`.
                """

        if (self.omega < 0) or (self.omega > 1):
            raise ValueError("Parameter omega must satify 0 <= omega =< 1.")

        start = time.time()
        self.M_nnz, self.M_rows, self.partition_table = fMPI.super_operator_extent(
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

        self.M_local_rows = self.partition_table[self.rank + 1] - self.partition_table[self.rank]
        self.M_n_row_starts = np.max([2, self.M_local_rows + 1])

        self.M_row_starts, self.M_col_indexes, self.M_values = fMPI.super_operator(
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
                self.M_nnz,
                self.M_n_row_starts,
                self.flock,
                self.MPI_communicator.py2f())

        finish = time.time()

        self.construction_time = finish - start

        start = time.time()
        self.M_num_rec_inds, self.M_rec_disps, self.M_num_send_inds, self.M_send_disps = fMPI.rec_a(
                self.M_rows,
                self.M_row_starts,
                self.M_col_indexes,
                self.partition_table,
                self.MPI_communicator.py2f())

        self.M_local_col_inds, self.M_rhs_send_inds = fMPI.rec_b(
                self.M_rows,
                np.sum(self.M_num_send_inds),
                self.M_row_starts,
                self.M_col_indexes,
                self.M_num_rec_inds,
                self.M_rec_disps,
                self.M_num_send_inds,
                self.M_send_disps,
                self.partition_table,
                self.MPI_communicator.py2f())
        finish = time.time()

        self.reconcile_time = finish - start

        start = time.time()

        self.one_norms, self.p = fMPI.one_norm_series(
                self.M_rows,
                self.M_row_starts,
                self.M_col_indexes,
                self.M_values,
                self.M_num_rec_inds,
                self.M_rec_disps,
                self.M_num_send_inds,
                self.M_send_disps,
                self.M_local_col_inds,
                self.M_rhs_send_inds,
                self.partition_table,
                self.MPI_communicator.py2f())
        finish = time.time()

        self.one_norm_time = finish - start

    def set_omega(self, omega):
        """Re-sets :math:`\omega` in the QSW master equation. This triggers :py:func:`construct_superoperator`.

        :param omega: :math:`\omega`
        :type omega: float
        """


        self.omega = omega
        self.construct_superoperator()

    def initial_state(self, state, name = None):
        """Sets the initial state of :math:`\\rho(0)`.

        :param state: One of the following:

            * (:math:`\\tilde{N}, \\tilde{N}`) or (N, N), complex,

            * :math:`\\tilde{N}` or N, complex,  Array describing the diagonal of :math:`\\rho(0)` with no inital inter-vertex coherences.

            * 'sources' - Evenly distribute the initial state over all specified sources.

            * 'even' - Evenly distribute the initial state over all sites, including specified sources and sinks.

        :param name: An identifying name for the initial state, for use when saving or loading states from a .qsw file.
        :type name: string, optional
        """

        self.initial_state_name = name

        if state is 'sources':
            self.rho = np.zeros((self.size[0]), dtype = np.complex128)
            self.rho[-(self.n_sinks + self.n_sources):-self.n_sinks] = np.float64(1)/np.float64(self.n_sources)
            self.rho = np.diag(self.rho)
            self.initial_state_name = 'sources'
            self.initial_state_gen = 'FreeQSW'

        elif state is 'even':
            self.rho = np.full((self.size[0]), np.float64(1)/np.float64(self.size[0]), dtype = np.complex128)
            self.rho = np.diag(self.rho)
            self.initial_state_name = 'even'
            self.initial_state_gen = 'FreeQSW'

        else:

            if not isinstance(state, np.ndarray):
                state = np.array(state, dtype = np.complex128)

            if len(state.shape) is 2:
                self.rho = np.concatenate((state, np.zeros((state.shape[0], self.pad))), axis = 1)
                self.rho = np.concatenate((self.rho, np.zeros((self.pad, self.rho.shape[1]))), axis = 0)
                self.initial_state_gen = 'User'

            elif len(state.shape) is 1:
                self.rho = np.diag(np.concatenate((state, np.zeros((self.pad,))), axis = 0))
                self.initial_state_gen = 'User'

        print(self.rho.shape)

        self.rho_v = fMPI.initial_state(
                self.M_local_rows,
                self.rho,
                self.rank,
                self.partition_table,
                self.MPI_communicator.py2f())

        if self.rank is 0:

            if self.filename is not None:

                output = h5py.File(self.filename + ".qsw", "a")

                if self.initial_state_name is None:
                    self.initial_state_name = 'rho ' + str(output['initial states'].attrs['counter'])
                    output['initial states'].attrs['counter'] += 1

                if self.initial_state_name not in output['initial states']:
                    output['initial states'].create_dataset(self.initial_state_name, data = self.rho, compression = "gzip")

    def File(self, filename, action = "a"):
        """Associate a .qsw file with the :class:`walk` object.

        If the file does not exist it is created, otherwise it is checked for consistency with the system of the walk. This file \
        may then be used to automatically save results obtained by :py:meth:`~freeqsw.MPI.walk.step` and :py:meth:`~freeqsw.MPI.walk.series`.

        :param filename: Filename of the .qsw file, or filename of the .qsw to be created. It must end with a .qsw extension.
        :type filename: string

        :param action: "a", open or create file. "w" create file, overwrite if it already exists.
        :type action: string
        """

        self.filename = str(filename)

        if self.rank is 0:
            io.set_file(
                    self.filename,
                    self.H,
                    self.L,
                    self.source_sites,
                    self.source_rates,
                    self.sink_sites,
                    self.sink_rates,
                    action)

    def step(self, t, target = None, save = False, name = None, precision = None):
        """Calculate :math:`\\rho(t)`.

        :param t: time.
        :type t: float

        :param target: MPI rank to receive :math:`\\rho(t)`.
        :type target: integer, optional

        :param save: If True, save the result to an associated .qsw file.
        :type save: boolean, optional

        :param name: Name under which to save the result, if None it defaults to 'step + <sequential number>'.
        :type name: string, optional

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision. Defaults to "dp".
        :type precision: string, optional

        :return rhot: If **target** is specified the :math:`\\rho(t)` is returned to the target rank.
        :rtype: (:math:`\\tilde{N}, \\tilde{N}`), complex

        .. warning::
            At minimum a **target** must be specified, or **save** must be 'True'.
        """

        if (target is None) and (save is False):
            raise ValueError("target= None and save=False. No target specified for series result")

        if precision is None:
            precision = "dp"

        start = time.time()

        rhot_v = fMPI.step(
                self.M_rows,
                self.M_row_starts,
                self.M_col_indexes,
                self.M_values,
                self.M_num_rec_inds,
                self.M_rec_disps,
                self.M_num_send_inds,
                self.M_send_disps,
                self.M_local_col_inds,
                self.M_rhs_send_inds,
                t,
                self.rho_v,
                self.partition_table,
                self.p,
                self.one_norms,
                self.MPI_communicator.py2f(),
                precision)

        end = time.time()

        self.step_time = end - start

        if save:

            if self.rank == 0:
                output_size = self.size
            else:
                output_size = 1

            rhot = fMPI.gather_step(
                    rhot_v,
                    self.partition_table,
                    0,
                    self.MPI_communicator.py2f(),
                    output_size)

            if self.rank is 0:

                output = h5py.File(self.filename + ".qsw", "a")

                if name is None:
                    name = 'step ' + str(output['steps'].attrs['counter'])
                    output['steps'].attrs['counter'] += 1

                step = output['steps'].create_dataset(name, data = rhot, compression = "gzip")
                step.attrs['omega'] = self.omega
                step.attrs['t'] = t
                step.attrs['initial state'] = self.initial_state_name
                step.attrs['initial state gen'] = self.initial_state_gen

                output.close()

        if target is not None:

            if (not save) or (target is not 0):

                if self.rank == 0:
                    output_size = self.size
                else:
                    output_size = 1

                rhot = fMPI.gather_step(
                        rhot_v,
                        self.partition_table,
                        target,
                        self.MPI_communicator.py2f(),
                        output_size)

            return rhot

    def series(self, t1, tq, steps, target = None, save = False, name = None, precision = None, chunk_size = None):
        """Calculate :math:`\\rho(t)` over :math:`[t_1, t_q]` at :math:`q` evenly spaced steps.

        :param t1: :math:`t_1`.
        :type t1: float

        :param tq: :math:`t_q`.
        :type tq: float

        :param steps: :math:`q`.
        :type steps: integer

        :param target: MPI rank to receive :math:`[\\rho(t_1),\\rho(t_2),...,\\rho(t_q)]`.
        :type target: integer, optional

        :param save: If True, save the result to an associated .qsw file.
        :type save: boolean, optional

        :param name: Name under which to save the result, if None it defaults to 'series + <sequential number>'.
        :type name: string, optional

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision. Defaults to "dp".
        :type precision: string, optional

        :return rhot: If **target** is specified :math:`[\\rho(t_1),\\rho(t_2),...,\\rho(t_q)]` is returned to the target rank.
        :rtype: (:math:`\\tilde{N}, \\tilde{N}`), complex, array

        .. warning::
            At minimum a **target** must be specified, or **save** must be 'True'.
        """

        if (target is None) and (save is False):
            raise ValueError("target= None and save=False. No target specified for series result")

        if precision is None:
            precision = "dp"

        if (chunk_size is not None) and (save is False):
            chunk_size = steps
            if self.rank == 0:
                warn('Chunked evaluation used only if save file is specified, defaulting to chunk_size=' + str(steps), UserWarning)

        if (chunk_size is None) and (save is False):
            chunk_size = steps
        elif (chunk_size is None) and (save is True):
            chunk_size = np.rint(float(40**9)/(16*self.H.shape[0]**2))

        chunks = int(np.ceil(steps/float(chunk_size)))

        h = np.float32(tq - t1)/np.float32(steps)
        chunk_steps = chunk_size

        for i in range(0,chunks):

            t1_chunk = i*chunk_size*h + t1

            if i == (chunks - 1):
                tq_chunk = tq
                chunk_steps = int(np.rint((tq_chunk - t1_chunk)/h))
            else:
                tq_chunk = (i + 1)*chunk_size*h + t1

            rhot_v_series = fMPI.series(
                    self.M_rows,
                    self.M_row_starts,
                    self.M_col_indexes,
                    self.M_values,
                    self.M_num_rec_inds,
                    self.M_rec_disps,
                    self.M_num_send_inds,
                    self.M_send_disps,
                    self.M_local_col_inds,
                    self.M_rhs_send_inds,
                    t1_chunk,
                    tq_chunk,
                    self.rho_v,
                    chunk_steps,
                    self.partition_table,
                    self.p,
                    self.one_norms,
                    self.MPI_communicator.py2f(),
                    precision)

            if save:

                if self.rank == 0:
                    output_size = self.size
                else:
                    output_size = 1

                rhot_series = fMPI.gather_series(
                        rhot_v_series,
                        self.partition_table,
                        0,
                        self.MPI_communicator.py2f(),
                        output_size)

                if (self.rank == 0) and (i == 0):

                    output = h5py.File(self.filename + ".qsw", "a")

                    if name is None:
                        name = 'series ' + str(output['series'].attrs['counter'])
                        output['series'].attrs['counter'] += 1

                    series = output['series'].create_dataset(
                            name,
                            data = rhot_series.T,
                            maxshape = (steps + 1, self.size[0], self.size[0]),
                            chunks = True,
                            compression = "gzip")

                    series.attrs['omega'] = self.omega
                    series.attrs['t1'] = t1
                    series.attrs['tq'] = tq
                    series.attrs['steps'] = steps
                    series.attrs['initial state'] = self.initial_state_name
                    series.attrs['initial state gen'] = self.initial_state_gen
                    output.close()

                elif self.rank == 0:

                    output = h5py.File(self.filename + ".qsw", "a")
                    series = output['series'][name]
                    series.resize(series.shape[0] + chunk_steps, axis = 0)
                    series[-chunk_steps:,:,:] = rhot_series.T[1:]

                    output.close()

        if target is not None:

            if (not save) or (target is not 0):

                if self.rank == 0:
                    output_size = self.size
                else:
                    output_size = 1

                rhot_series = fMPI.gather_series(
                        rhot_v_series,
                        self.partition_table,
                        target,
                        self.MPI_communicator.py2f(),
                        output_size)

                return rhot_series.T

            else:

                if self.rank == target:
                    output = h5py.File(self.filename + ".qsw", "r")
                    rhot_series = output['series'][name][:]
                    output.close()

                    return rhot_series

def load_system(filename, MPI_communicator):
    """
    Load :math:`H`, :math:`L`, sources and sinks stored in a .qsw file and distribute this over an MPI communicator.

    :param filename: Name of the .qsw file.
    :type filename: string

    :param MPI_communicator: MPI communicator.

    :return: :math:`H`, :math:`L`, sources and sinks.
    :rtype: (N, N), complex; ([M], [M]), (integer, float), tuple; ([M], [M]), (integer, float), tuple.
    """

    rank = MPI_communicator.Get_rank()

    if rank == 0:
        system = io.File(str(io._qsw_extension(filename)))
        H = system.H()
        L = system.L()
        sources = system.sources()
        sinks = system.sinks()
        size = H.shape[0]
        system.close()
    else:
        size = None
        system = None
        sources = (None, None)
        sinks = (None, None)

    size = MPI_communicator.bcast(size, root = 0)
    sources = MPI_communicator.bcast(sources, root=0)
    sinks = MPI_communicator.bcast(sinks, root=0)

    if rank != 0:
        H = sp.csr_matrix((size,size), dtype = np.complex128)
        L = sp.csr_matrix((size,size), dtype = np.complex128)

    H.data = MPI_communicator.bcast(H.data, root = 0)
    H.indices = MPI_communicator.bcast(H.indices, root = 0)
    H.indptr = MPI_communicator.bcast(H.indptr, root = 0)

    L.data = MPI_communicator.bcast(L.data, root = 0)
    L.indices = MPI_communicator.bcast(L.indices, root = 0)
    L.indptr = MPI_communicator.bcast(L.indptr, root = 0)

    return H, L, sources, sinks
