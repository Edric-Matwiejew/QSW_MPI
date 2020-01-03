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
import qsw_mpi.foperators as foperators

def check_indices(G):
    if not G.has_sorted_indices:
        G.sort_indices()

def transition(gamma, G):
    """Generate the graph transition matrix from an adjacency matrix.

    :param gamma: :math:`\gamma`.
    :type gamma: float

    :param G: :math:`G`.
    :type G: SciPy complex CSR matrix

    :return M: :math:`M`
    :rtype: SciPy complex CSR matrix
    """

    if gamma < 0:
        raise ValueError("0 < gamma")

    check_indices(G)

    M_indptr, M_indices, M_data, M_nnz = foperators.graph(gamma, G.indptr, G.indices, G.data)

    M = sp.csr_matrix((M_data[0:M_nnz], M_indices[0:M_nnz], M_indptr), G.shape, dtype = np.complex128)

    return M

def lindblad(G):
    """Generate a Lindblad operator matrix from a graph transition matrix.

    :param G: :math:`G`.
    :type G: SciPy complex CSR matrix

    :return L: :math:`L`
    :rtype: SciPy complex CSR matrix
    """
    check_indices(G)

    L_indptr, L_indices, L_data, L_nnz = foperators.site_lindblads(G.indptr, G.indices, G.data)

    L = sp.csr_matrix((L_data[0:L_nnz], L_indices[0:L_nnz], L_indptr), G.shape)

    return L

def symmetrise(G):
    """Symmetrise a CSR matrix with a symmetric structure.

    G is altered so as to satisfy :math:`G_{ij} = max(G_{ij}, G_{ji})`.

    :param G: CSR matrix to symmetrise.
    :type G: SciPy complex CSR matrix.

    .. warning::
        The matrix must be structurally symmetric in terms of its non-zero entires.
    """
    check_indices(G)

    foperators.symmetrise(G.indptr, G.indices, G.data)
