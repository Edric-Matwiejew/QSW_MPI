import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import h5py
import freeqsw.operators as operators
import freeqsw.fMPI as fMPI
import freeqsw.io as io

def load_walk(
        omega,
        filename,
        MPI_communicator):
    """ Initialize a previously created :class:`walk` using a .qsw file.

    :param omega: Interpolation paramater (:math:`\omega`) \
            in the stochastic walk master equation, \
            must satisfy :math:`0 \geq \omega \leq 1`.
    :type omega: numpy.real64

    :param filename: Target file. Must have a .qsw extention and not be currently open.
    :type filename: str

    :param MPI_communicator: An MPI communicator created using mpi4py.
    :type MPI_communicator: MPI intracommunicator

    :return: A :class:`walk` object initialized using the input systems's H, L, \
            sources, and sinks data. The :class:`walk.file` parameter is set to the \
            input file.
    """
    H, L, sources, sinks = io.load_system(filename, MPI_communicator)

    walk = walk(omega, H, L, MPI_communicator, sources, sinks)
    walk.file(filename, action = "a")

    return walk

class walk(object):
    """The :class:`walk` object handels the creation of distributed superoperators, density matrices and system propagation in an MPI environment. It also enables the saving of results and input data to a .qsw file. This object and its contianing methods will only function correctly if called within an MPI instance, as shown in the examples.

    :param omega: Interpolation paramater (:math:`\omega`) \
        in the stochastic walk master equation, \
        must satisfy :math:`0 \geq \omega \leq 1`.
    :type omega: numpy.real64

    :param H: The system Hamiltonian, a Hermitian matrix.
    :type H: numpy.complex128

    :param L: Matrix resulting from summation of the system Lindblad operators.\
            May be non-hermitian.
    :type L: numpy.complex128

    :param MPI_communicator: An MPI communicator created using mpi4py.
    :type MPI_communicator: MPI intracommunicator

    :param sources: Tuple containing a list of  absorption sites and the \
            corresponding absorption rates.
    :type sources: (numpy.int32, numpy.float64), optional

    :param sinks: Tuple containing a list of the emission sites and the\
            corresponding absorption rates.
    :type sinks: (numpy.int32, numpy.float64), optional

    :return: A distributed :class:`walk` object containing an initialized \
            superoperator.


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
            self.source_sites = sources[0][source_idx]
            self.source_rates = sources[1][source_idx]

            sink_idx = np.argsort(sinks[0])
            self.sink_sites = sinks[0][sink_idx]
            self.sink_rates = sinks[1][sink_idx]

        elif self.source_sites is not None:

            self.pad = len(self.source_sites)
            self.n_sources = len(self.source_sites)
            self.n_sinks = 0
            self.sink_sites = [-1]
            self.sink_rates = [0]

            source_idx = np.argsort(sources[0])
            self.source_sites = sources[0][source_idx]
            self.source_rates = sources[1][source_idx]

            sink_idx = np.argsort(sinks[0])
            self.sink_sites = sinks[0][sink_idx]
            self.sink_rates = sinks[1][sink_idx]

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

        self.construct_superoperator()

    def construct_superoperator(self):
        """A distributed superoperator is constructed on :class:`walk` initialization\
                or on calling :py:meth:`~freeqsw.MPI.walk.set_omega`. One-norms \
                of the superoperator power series are calculated for using in \
                deriving the optimal matrix exponentiation parameters when calling \
                :py:meth:`~freeqsw.MPI.walk.step` and :py:meth:`~freeqsw.MPI.walk.series`.
                """

        if (self.omega < 0) or (self.omega > 1):
            raise ValueError("Parameter omega must satify 0 <= omega =< 1.")

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
                self.M_local_rows,
                self.flock,
                self.MPI_communicator.py2f())

        self.M_num_rec_inds, self.M_rec_disps, self.M_num_send_inds, self.M_send_disps = fMPI.rec_a(
                self.M_rows,
                self.M_row_starts,
                self.M_col_indexes,
                self.partition_table,
                self.MPI_communicator.py2f())

        self.M_local_col_inds, self.M_rhs_send_inds = fMPI.rec_b(
                self.M_rows,
                np.sum(self.M_num_rec_inds),
                self.M_row_starts,
                self.M_col_indexes,
                self.M_num_rec_inds,
                self.M_rec_disps,
                self.M_num_send_inds,
                self.M_send_disps,
                self.partition_table,
                self.MPI_communicator.py2f())

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

    def set_omega(self, omega):
        """Re-sets the interpolation paramater (:math:`\omega`) \
            in the stochastic walk master equation. This triggers \
            :py:func:`construct_superoperator`.

        :param omega: Must satisfy :math:`0 \geq \omega \leq 1`.
        :type omega: numpy.real64
        """


        self.omega = omega
        self.construct_superoperator()

    def initial_state(self, state, name = None):
        """Sets the initial state of the density operator (:math:`\\rho`).

        :param state:
        :type state: np.complex128 or str

        Possible input are:

        * A 2D-array with shape equal to H.shape or shape equal to H.shape \
                + the number of sources and sinks.

        * A 1D-array with length equal to the diagonal of H or length equal to the \
                diagonal of H + the number of sources and sinks. This describes an \
                inital state with no coherences.

        * 'sources' - evenly distribute the initial state over \
                all specified sources.

        * 'even' - evenly distribute the initial state over all \
                sites, including an specified sources and sinks.

        If needed the array will be padded to match the shape of the system operators \
                after the addtion of absorption and emission sites (sources and sinks).

        :param name: An identifying name for the initial state, for use when saving \
                or loading states from a .qsw file.
        :type name: str, optional
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
                self.rho = np.concatenate((state, np.zeros((self.pad, state.shape[0]))), axis = 0)
                self.initial_state_gen = 'User'

            elif len(state.shape) is 1:
                self.rho = np.diag(np.concatenate((state, np.zeros((self.pad,))), axis = 0))
                self.initial_state_gen = 'User'

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

        If the file does not exist it is created, otherwise it is \
        checked for consistency with the system of the walk. This file \
        may then be used to automatically save results obtained by \
        :py:meth:`~freeqsw.MPI.walk.step` and :py:meth:`~freeqsw.MPI.walk.series`.

        :param filename: Filename of the .qsw file, or filename of the .qsw \
                to be created. It must end with a .qsw extension.
        :type filename: str

        :param action: "a", open or create file. "w" create file, overwrite if it already exists.
        :type action: str
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
        """Calculate the state of the system at time t.

        :param t: The time, must statisfy :math:`0 \ge t`.
        :type t: numpy.float64

        :param target: MPI rank to receive the calculated density operator, :math:`\\rho(t)`.
        :type target: int, optional

        :param save: If True, save the result to an associated .qsw file.
        :type save: True or False, optional

        :param name: Name under which to save the result, if None it defaults to 'step + <sequential number>'.
        :type name: str, optional

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision. Defaults to "dp".
        :type precision: str, optional

        :return rhot: If **target** is specified the :math:`\\rho(t)` is returned to the target rank.
        :rtype: numpy.complex128

        .. warning::
            At minimum a **target** must be specified, or **save** must be 'True'.
        """

        if precision is None:
            precision = "dp"


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

        if save:

            rhot = fMPI.gather_step( rhot_v, \
                                    self.partition_table, \
                                    0, \
                                    self.MPI_communicator.py2f(), \
                                    self.size)

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

                rhot = fMPI.gather_step( rhot_v, \
                                        self.partition_table, \
                                        target, \
                                        self.MPI_communicator.py2f(), \
                                        self.size)
            return rhot

    def series(self, t1, t2, steps, target = None, save = False, name = None, precision = None):
        """Calculate the state of the system over the time interval t1 to t2.

        :param t1: The starting time, must statisfy :math:`0 \ge t`.
        :type t1: numpy.float64

        :param t2: The ending time, must statisfy :math:`t1 \ge t2`.
        :type t2: numpy.float64

        :param steps: The number of evenly spaced points in time to calculate.

        :param target: MPI rank to receive the calculated density operator, :math:`\\rho(t)`.
        :type target: int, optional

        :param save: If True, save the result to an associated .qsw file.
        :type save: True or False, optional

        :param name: Name under which to save the result, if None it defaults to 'series + <sequential number>'.
        :type name: str, optional

        :param precision: Target precision, accepts "sp" for single precision and "dp" for double precision. Defaults to "dp".
        :type precision: str, optional

        :return rhot: If **target** is specified the :math:`\\rho(t)` is returned to the target rank.
        :rtype: numpy.complex128

        .. warning::
            At minimum a **target** must be specified, or **save** must be 'True'.
        """


        if precision is None:
            precision = "dp"

        rhot_v_series = fMPI.series( self.M_rows, \
                                    self.M_row_starts, \
                                    self.M_col_indexes, \
                                    self.M_values, \
                                    self.M_num_rec_inds, \
                                    self.M_rec_disps, \
                                    self.M_num_send_inds, \
                                    self.M_send_disps, \
                                    self.M_local_col_inds, \
                                    self.M_rhs_send_inds, \
                                    t1, \
                                    t2, \
                                    self.rho_v, \
                                    steps, \
                                    self.partition_table, \
                                    self.p, \
                                    self.one_norms, \
                                    self.MPI_communicator.py2f(), \
                                    precision)

        if save:

            rhot_series = fMPI.gather_series(rhot_v_series, \
                                            self.partition_table, \
                                            0, \
                                            self.MPI_communicator.py2f(), \
                                            self.size)

            if self.rank is 0:

                output = h5py.File(self.filename + ".qsw", "a")

                if name is None:
                    name = 'series ' + str(output['series'].attrs['counter'])
                    output['series'].attrs['counter'] += 1

                series = output['series'].create_dataset(name, data = rhot_series.T, compression = "gzip")
                series.attrs['omega'] = self.omega
                series.attrs['t1'] = t1
                series.attrs['t2'] = t2
                series.attrs['steps'] = steps
                series.attrs['initial state'] = self.initial_state_name
                series.attrs['initial state gen'] = self.initial_state_gen

                output.close()

        if target is not None:

            if (not save) or (target is not 0):

                rhot_series = fMPI.gather_series(rhot_v_series, \
                                                self.partition_table, \
                                                target, \
                                                self.MPI_communicator.py2f(), \
                                                self.size)

            return rhot_series.T

