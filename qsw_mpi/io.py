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
import h5py
from os import path

"""
Non-paralel input/output operations.
"""

def _qsw_extension(filename):
    """
    Ensures input and output file names end with a .qsw extension.
    """
    if path.splitext(filename)[-1] == '.qsw':
        return filename
    else:
        return filename + '.qsw'

class File(object):

    """ :class:`File` handels the creation of .qsw files and the accessing of data contained within existing .qsw files. It also provides convenience functions to list the :meth:`~qsw_mpi.MPI.walk.step` and :meth:`~qsw_mpi.MPI.series` data contained within an existing .qsw file.

    :param filename: File name of the .qsw file to open or create.
    :type filename: string
    """

    def __init__(self, filename):
        self.File = h5py.File(_qsw_extension(filename),"a")

    def load_csr(self, group):
        """
        Retrive an (:math:`\\tilde{N}, \\tilde{N}`) SciPy CSR sparse matrix from a .qsw file group.

        :param group: Name of the target group in the .qsw file.
        :type group: string

        :returns A:
        :rtype: SciPy complex CSR matrix.
        """
        csr = sp.csr_matrix(
                (np.array(self.File[group]['data'], dtype = np.complex128),
                np.array(self.File[group]['indices'], dtype =np.int32),
                np.array(self.File[group]['indptr'], dtype =np.int32)),
                shape = (self.File[group].attrs['size'], self.File[group].attrs['size']))

        return csr

    def H(self):
        """
        Returns the graph Hamiltonian associated with the .qsw file.

        :return: H
        :type: SciPy complex CSR matrix
        """
        return self.load_csr('H')

    def L(self):
        """
        Returns the Lindblad operator matrix associated with the .qsw file.

        :return: H
        :type: Scipy complex CSR matrix
        """
        return self.load_csr('L')

    def sources(self):
        """
        Returns sources associated with the .qsw file.

        :return:  A tuple containing two arrays. `source_sites` contains the graph node to which each source is connected and `source_rates` contains the transition rate from each source.
        :rtype: ([M], [M]), (integer, float), tuple
        """
        if self.File.attrs['# sources'] > 0:
            return (np.array(self.File['sources']['source sites'], dtype = np.int32), \
                    np.array(self.File['sources']['source rates'], dtype = np.float64))
        else:
            return (None, None)

    def sinks(self):
        """
        Returns sinks associated with the .qsw file.

        :return: A tuple containing two arrays. `sink_sites` contains the graph node to which each source is connected and `sink_rates` contains the transition rate to each sink.
        :rtype: ([M], [M]), (integer, float), tuple
        """

        if self.File.attrs['# sinks'] > 0:
            return (np.array(self.File['sinks']['sink sites'], dtype = np.int32), \
                    np.array(self.File['sinks']['sink rates'], dtype = np.float64))
        else:
            return (None, None)

    def initial_state(self, state):
        """
        Returns the initial state of the density operator, :math:`\\rho(0)`, associated with a .qsw file.


        :param state: Name of the initial state.
        :type state: string

        :return: :math:`\\rho(0)`.
        :rtype: (:math:`\\tilde{N}, tilde{N}`), complex
        """

        if self.File['initial states gen'][str(name)] is 'FreeQSW':
            raise NameError("FreeQSW generated initial states are not saved to file.")

        return np.array(self.File['initial states'][str(name)], dtype = np.complex128)

    def step(self, name):
        """
        Returns a saved :math:`\\rho(t)`.

        :param name: Name of the step.
        :type name: string

        :return: :math:`\\rho(t)`.
        :rtype: (:math:`\\tilde{N}, \\tilde{N}`)
        """
        return np.array(self.File['steps'][str(name)], dtype = np.complex128)

    def series(self, name):
        """
        Returns a saved time series, :math:`[\\rho(t_1),\\rho(t_2), ..., \\rho(t_q)]`.

        :param name: Name of the series.
        :type name: string

        :return: :math:`[\\rho(t_1),\\rho(t_2),..., \\rho(t_q)]`
        :rtype: (:math:`\\tilde{N}, \\tilde{N}`), complex, array
        """
        return np.array(self.File['series'][str(name)], dtype = np.complex128)

    def list_steps(self, spacing = 15):
        """
        Gives a summary of the :math:`\\rho(t)` stored in the .qsw file. `list_steps` prints the step names, initial state names (:math:`\\rho(t_0)`), omega (:math:`\\omega`) and the simulation times (:math:`t`) in a human-readable format.

        :param spacing: Optional, default = 15. The width of each column in number of characters.
        :type spacing: integer
        """
        titles = ["step name", "initial state", "omega", "t"]

        step_names = []
        inital_states = []
        omegas = []
        t = []

        for item in self.File['steps']:
            step_names.append(str(item))
            inital_states.append(str(self.File['steps'][str(item)].attrs['initial state']))
            omegas.append(float(self.File['steps'][str(item)].attrs['omega']))
            t.append(float(self.File['steps'][str(item)].attrs['t']))

        data = [titles] + list(zip(step_names, inital_states, omegas, t))

        for i, d in enumerate(data):
            line = '|'.join(str(x).ljust(spacing) for x in d)
            print(line)
            if i == 0:
                print('-' * len(line))

    def list_series(self, spacing = 15):
        """
        Gives a summary of the :math:`\\rho_{series} = [\\rho(t_1), ..., \\rho(t_q)]` stored in the .qsw file. `list_series` prints the series names, initial state names (:math:`\\rho(t_1)`), omega (:math:`\\omega`) and the simulation times (:math:`t_1, t_q`) in a human-readable format.

        :param spacing: Optional, default = 15. The width of each column in number of characters.
        :type spacing: integer
        """

        titles = ["series name", "initial state", "omega", "steps", "t1", "tq"]

        series_names = []
        inital_states = []
        omegas = []
        steps = []
        t1 = []
        tq = []

        for item in self.File['series']:
            series_names.append(str(item))
            inital_states.append(str(self.File['series'][str(item)].attrs['initial state']))
            omegas.append(float(self.File['series'][str(item)].attrs['omega']))
            steps.append(int(self.File['series'][str(item)].attrs['steps']))
            t1.append(float(self.File['series'][str(item)].attrs['t1']))
            tq.append(float(self.File['series'][str(item)].attrs['tq']))

        data = [titles] + list(zip(series_names, inital_states, omegas, steps, t1, tq))

        for i, d in enumerate(data):
            line = '|'.join(str(x).ljust(spacing) for x in d)
            print(line)
            if i == 0:
                print('-' * len(line))

    def system_attributes(self):
        """
        Return meta-data contained with a .qsw file.

        :return: {"size": :math:`\\tilde{N}` (integer), \
                "H nnz": non-zeros in :math:`H` (integer), \
                "L nnz": non-zeros in :math:`L` (integer), \
                "# sources": number of sources (integer), \
                "# sinks": number of sinks (integer)}
        :rtype: dictionary
        """
        system = {
                "size":int(self.File['H'].attrs['size']),
                "H nnz":int(self.File['H'].attrs['nnz']),
                "L nnz":int(self.File['L'].attrs['nnz']),
                "# sources":int(self.File.attrs["# sources"]),
                "# sinks":int(self.File.attrs["# sinks"])}
        return system

    def close(self):
        """
        Close the opened .qsw file.
        """
        self.File.close()

def set_file(
        filename,
        H,
        L,
        source_sites,
        source_rates,
        sink_sites,
        sink_rates,
        action = "a"):

    """
    Creates a .qsw file on disk storing the information needed to recreate a given system. If the file already exists, it can be set as the destination for additional simulations.

    :param filename: Name of the .qsw file.
    :type filename: string

    :param H: :math:`H`
    :type H: (N, N), complex

    :param L: :math:`L`
    :type L: (N, N), complex

    :param source_sites: Vertex in :math:`G` to which a source is attached.
    :type source_sites: integer

    :param source_rates: Transition rates of attached sources.
    :type source_rates: float

    :param sink_sites: Vertex in :math:`G` to which a sink is attached.
    :type sink_sites: integer

    :param sink_rates: Transition rates of attached sinks.
    :type sink_rates: float

    :param action: Optional, default = "a". h5py file creation flag.
    :type action: string
    """
    output = h5py.File(_qsw_extension(filename), action)

    if 'H' not in output:
        H_out = output.create_group('H')
        H_out.attrs['size'] = H.shape[0]
        H_out.attrs['nnz'] = H.count_nonzero()
        H_out.create_dataset('data', data = H.data)
        H_out.create_dataset('indptr', data = H.indptr)
        H_out.create_dataset('indices', data = H.indices)

    elif (int(H.shape[0])!= int(output['H'].attrs['size'])) \
            or (int(H.count_nonzero()) != int(output['H'].attrs['nnz'])):
        raise NameError('Stored H does not match with the current walk.')

    if 'L' not in output:

        L_out = output.create_group('L')
        L_out.attrs['size'] = L.shape[0]
        L_out.attrs['nnz'] = L.count_nonzero()
        L_out.create_dataset('data', data = L.data)
        L_out.create_dataset('indptr', data = L.indptr)
        L_out.create_dataset('indices', data = L.indices)

    elif (int(L.shape[0]) != int(output['L'].attrs['size'])) \
            or (int(L.count_nonzero()) != int(output['L'].attrs['nnz'])):
        raise NameError('Stored L does not match with the current walk.')

    if 'sources' not in output:

        if source_sites[0] is not -1:
            output.attrs['sources'] = True
            output.attrs['# sources'] = source_sites.shape[0]
            sources = output.create_group('sources')
            sources.create_dataset('source sites', data = source_sites)
            sources.create_dataset('source rates', data = source_rates)
        else:
            output.attrs['# sources'] = 0

    elif source_sites.shape[0] != int(output.attrs['# sources']):
            raise NameError('Stored source sites do not match with the current walk.')

    if 'sinks' not in output:
        if sink_sites[0] is not -1:
            output.attrs['sinks'] = True
            output.attrs['# sinks'] = sink_sites.shape[0]
            sinks = output.create_group('sinks')
            sinks.attrs['sink sites'] = len(sink_sites)
            sinks.create_dataset('sink sites', data = sink_sites)
            sinks.create_dataset('sink rates', data = sink_rates)
        else:
            output.attrs['# sinks'] = 0


    elif source_sites.shape[0] != int(output.attrs['# sinks']):
            raise NameError('Stored sink sites do not match with the current walk.')

    if 'initial states' not in output:
        output.create_group('initial states')
        output['initial states'].attrs['counter'] = 1

    if 'steps' not in output:
        output.create_group('steps')
        output['steps'].attrs['counter'] = 1

    if 'series' not in output:
        output.create_group('series')
        output['series'].attrs['counter'] = 1

    output.close()

