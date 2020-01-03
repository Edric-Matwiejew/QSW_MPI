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

"""
Non-parallel extraction and manpulation of data from QSW simulations.
"""

def populations(rho = None, File = None, step_name = None, series_name = None):
    """
    Return the populations, :math:`\\text{trace}(\\rho(t))`, of :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]`.

    This can be done in two ways.

        #. Pass :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]` as a numpy array via the `rho` parameter.

        #. Access :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]` from a .qsw file by passing a :class:`File` object by the `File` parameter and its name by the `step_name` or `series_name` parameter.

    :param rho: :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]`
    :type rho: (:math:`\\tilde{N}, \\tilde{N}`), complex, array

    :param File: .qsw filename.
    :type File: string

    :param step_name: Name of a step.
    :type step_name: string

    :param series_names: Name of a series.
    :type series_names: string

    :return: :math:`\\text{trace}(\\rho(t))` or :math:`[\\text{trace}(\\rho(t_1)), ..., \\text{trace}(\\rho(t_q))]`
    :rtype: :math:`\\tilde{N}`, float, array
    """
    if rho is not None:

        if rho.ndim is 2:
            site_pops = np.real(np.diag(rho))
            return site_pops

        elif rho.ndim is 3:
            site_pops = np.array([np.real(np.diag(rho_step)) for rho_step in rho])
            return site_pops

    if File is not None:

        if step_name is not None:
            site_pops = np.real(np.diag(np.array(File.File['steps'][str(step_name)], dtype = np.complex128)))
            return site_pops

        elif series_name is not None:
            steps = File.File['series'][str(series_name)].attrs['steps'] + 1
            site_pops = np.array(
                    [np.real(np.diag(np.array(File.File['series'][str(series_name)][i,:,:],
                    dtype = np.complex128))) for i in range(steps)])
            return site_pops

def coherences(rho = None, File = None, step_name = None, series_name = None):
    """
    Return the inter-vertex coherences, :math:`\\\\rho(t)_{ij}, i > j`, of a :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]`.

    This can be done in two ways.

        #. Pass :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]` as a numpy array via the `rho` parameter.

        #. Access :math:`\\rho(t)` or :math:`[\\rho(t_1), \\rho(t_2), ..., \\rho(t_q)]` from a .qsw file by passing a :class:`File` object by the `File` parameter and its name by the `step_name` or `series_name` parameter.

    :param rho: :math:`\\rho(t)` or :math:`[\\rho(t_0), ..., \\rho(t_q)]`
    :type rho: (:math:`\\tilde{N}, \\tilde{N}`), complex, array

    :param File: .qsw filename.
    :type File: string

    :param step_name: Name of a step.
    :type step_name: string

    :param series_names: Name of a series.
    :type series_names: string

    :return: * :math:`[(i,j), ..., (\\tilde{N}, \\tilde{N} - 1)]`
             * :math:`\\rho(t))_{ij}, i>j` or :math:`[\\rho(t_0)_{ij}, ..., \\rho(t_q)_{ij}], i>j`
    :rtype: :math:`\\tilde{N}(\\tilde{N} - 1)/2`, complex, array
    """


    if rho is not None:

        if rho.ndim is 2:

            node_pairs = np.triu_indices_from(rho, k=1)
            cohs = np.asarray(abs(rho[node_pairs]))

            return node_pairs, cohs

        elif rho.ndim is 3:

            node_pairs = np.triu_indices_from(rho[0], k=1)
            cohs = np.asarray(abs(rho[0][node_pairs]))
            it = iter(rho)
            next(it,None)

            for step in it:

                cohs = np.vstack((cohs,np.asarray(abs(step[node_pairs]))))

            return node_pairs, cohs

    if File is not None:

        if step_name is not None:

            rho = np.array(File.File['steps'][str(step_name)], dtype = np.complex128)
            node_pairs = np.triu_indices_from(rho, k=1)
            cohs = np.asarray(abs(rho[node_pairs]))

            return node_pairs, cohs

        elif series_name is not None:

            steps = File.File['series'][str(series_name)].attrs['steps'] + 1
            rho = np.array(File.File['series'][str(series_name)][0], dtype = np.complex128)
            node_pairs = np.triu_indices_from(rho, k=1)
            cohs = np.asarray(abs(rho[node_pairs]))

            it = iter(range(steps))
            next(it, None)

            for i in it:
                cohs = np.vstack(
                        (cohs,np.asarray(abs(np.array(File.File['series'][str(series_name)][i,:,:],
                            dtype = np.complex128)[node_pairs]))))
            return node_pairs, cohs
