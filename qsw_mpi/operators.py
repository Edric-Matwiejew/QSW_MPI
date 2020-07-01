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
import scipy as sp
from scipy.special import factorial
import qsw_mpi.foperators as foperators
import copy

I = np.complex(0,1)

def check_indices(G):
    """
    Ensure that as SciPY CSR matrix has column ordered `indices` and `data` arrays.

    :param G: SciPy CSR matrix
    """
    if not G.has_sorted_indices:
        G.sort_indices()

def import_superoperator(File, operatorname):
    """
    Import a :math:`\\tilde{\mathcal{L}}` from a HDF5 file.

    .. Note::
        This method does not distrubted :math:`\\tilde{\mathcal{L}}`.

    :param File: H5Py :class:`File` object.

    :param operatorname: Group containing the :math:`\\tilde{\mathcal{L}}` CSR arrays.
    :type operatorname: string

    :returns: Superoperator, :math:`\\tilde{\mathcal{L}}`.
    :rtype: :math:`\\tilde{N}^2 \\times \\tilde{N}^2` SciPy CSR matrix
    """
    M = sp.sparse.csr_matrix(
            (
                np.array(File[operatorname + '/data']),
                np.array(File[operatorname + '/indices']),
                np.array(File[operatorname + '/indptr'])),
                File[operatorname].attrs['dimensions']
            )
    return M

def transition(gamma, G):
    """
    Generate the graph transition matrix (see Equation :eq:`eq:generator_matrix`) from a CSR adjacency matrix.

    :param gamma: Transition rate, :math:`\gamma`.
    :type gamma: float

    :param G: Adjacency matrix, :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix

    :return: The transition matrix of :math:`G`, :math:`M`
    :rtype: :math:`N \\times N` SciPy CSR matrix
    """

    if gamma < 0:
        raise ValueError("0 < gamma")

    check_indices(G)

    M_indptr, M_indices, M_data, M_nnz = foperators.graph(gamma, G.indptr, G.indices, G.data)

    M = sp.sparse.csr_matrix((M_data[0:M_nnz], M_indices[0:M_nnz], M_indptr), G.shape, dtype = np.complex128)

    return M

def symmetrise(G):
    """Symmetrise a CSR matrix with a symmetric structure.

    G is altered so as to satisfy :math:`G_{ij} = max(G_{ij}, G_{ji})`.

    :param G: CSR matrix to symmetrise.
    :type G: SciPy real CSR matrix.

    .. warning::
        The matrix must be structurally symmetric in terms of its non-zero entires.
    """

    check_indices(G)
    G_sym = copy.deepcopy(G)

    foperators.symmetrise(G_sym.indptr, G_sym.indices, G_sym.data)

    return G_sym

### L-QSW Methods ###

def local_lindblads(G):
    """
    Generate a Lindblad operator matrix from a graph transition matrix or adjacency matrix for use with the :class:`~qsw_mpi.MPI.LQSW` class (see Equation :eq:`eq:condensed_lindblads`).

    :param G: Adjacency matrix :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix

    :return L: :math:`M_L`.
    :rtype: :math:`N \\times N` SciPy CSR matrix
    """
    check_indices(G)

    L_indptr, L_indices, L_data, L_nnz = foperators.site_lindblads(G.indptr, G.indices, G.data)

    L = sp.sparse.csr_matrix((L_data[0:L_nnz], L_indices[0:L_nnz], L_indptr), G.shape)

    return L

def markov_chain(G):
    """
    Generate the cannonical Markov chain transition matrix (see Equation :eq:`eq:markov_chain`) of a CSR adjacency matrix for use with the :class:`~qsw_mpi.MPI.LQSW` class (see Equation :eq:`eq:condensed_lindblads`).


    :param G: Adjacency matrix, :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix.

    :return: The cannonical Markov chain transition matrix of :math:`G`, :math:`C`.
    :rtype: :math:`N \\times N` SciPy CSR matrix
    """

    outdegrees = np.zeros(G.shape[0],dtype=np.complex128)

    for i in range(G.shape[0]):
        for j in range(G.indptr[i], G.indptr[i+1]):
            outdegrees[i] += G.data[j]

    row_ind = []
    col_ind = []
    data = []
    for i in range(G.shape[0]):
        for j in range(G.shape[0]):
            if G[i,j] != 0:
                row_ind.append(i)
                col_ind.append(j)
                data.append((1.0/outdegrees[i]))

    return sp.sparse.csr_matrix((data,(row_ind, col_ind)), shape = G.shape)

### G-QSW Demoralization Methods ###

def __perm(n, i):
    """
    Generate the :math:`i^\\text{th}` permutation. Used to create column permutation of the Fourier matrix.
    """
    fact = factorial([i for i in range(n)], exact = True)

    perm = []
    for k in range(n):
        perm.append(i // fact[n - k - 1])
        i = i % fact[n - k - 1]

    for k in range(n - 1, 0, -1):
        for j in range(k - 1, -1, -1):
            if perm[j] <= perm[k]:
                perm[k] += 1

    return perm

def __fourier_matrix(n, perm_i = 0):

    """
    Returns a Fourier matrix of size :math:`n` with columns in the :math:`i^\\text{th}` permutation.
    """
    if ((perm_i == 0) or (n == 1)):
        col_inds = [i for i in range(n)]
    else:
        col_inds = __perm(n, perm_i)

    f = np.zeros((n,n), dtype = np.complex128)

    for i in range(n):
        for j in range(n):
            f[j, col_inds[i]] = np.exp((2 * np.pi * I * i * j) / np.float64(n))

    return f

def nm_vsets(G):
    """
    Creates vertex subspaces, :math:`V^D`, as per :ref:`Step 1 <demoral>` of the graph demoralisation procedure.

    :param G: Adjacency matrix, :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix

    :returns: Vertex subspaces, :math:`V^D`.
    :rtype: integer arrays
    """
    indegrees = []
    for i in range(G.shape[0]):
        indegrees.append(G.indptr[i+1]-G.indptr[i])

    vsets = [[] for i in range(G.shape[0])]

    index = 0
    for i in range(G.shape[0]):
        if indegrees[i] > 0:
            for j in range(indegrees[i]):
                vsets[i].append(index)
                index += 1
        else:
            vsets[i].append(index)
            index += 1

    return vsets

def import_vsets(File, vsetsname):
    """
    Import :math:`V^D` (see :ref:`Step 1 <demoral>` of the demoralisation procedure) from a HDF5 file.

    :param File: H5Py :class:`File` object.

    :returns: Vertex subspaces, :math:`V^D`.
    :rtype: integer arrays
    """

    vsets = []
    for key in File[vsetsname].keys():
        vsets.append(np.array(File[vsetsname][str(key)]))
    return vsets

def nm_G(G, vsets):
    """
    Returns the demoralised graph :math:`G^D`, as per :ref:`step 2 <demoral>` of the graph demoralisation procedure.

    :param G: Adjacency matrix, :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :returns: The demoralised graph, :math:`G^D`.
    :rtype: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR matrix
    """
    row_ind = []
    col_ind = []
    data = []

    sub_counts = np.zeros((len(vsets),len(vsets)), dtype = np.complex128)
    for i in range(G.shape[0]):
        for j in range(G.shape[0]):
            if G[i,j] != 0:
                for k in range(len(vsets[i])):
                    for l in range(len(vsets[j])):
                        sub_counts[i,j] += G[i,j]

    for i in range(G.shape[0]):
        for j in range(G.shape[0]):
            if G[i,j] != 0:
                for k in range(len(vsets[i])):
                    for l in range(len(vsets[j])):
                        row_ind.append(vsets[i][k])
                        col_ind.append(vsets[j][l])
                        data.append((1.0/np.sqrt(sub_counts[i,j])))

    shape = (vsets[-1][-1] + 1, vsets[-1][-1] + 1)
    return sp.sparse.csr_matrix((data,(row_ind, col_ind)), shape = shape)

def nm_L(G, vsets, perm_i = 0):
    """
    Returns a global Lindblad operator, :math:`L^D`, capable of generating a non-moralising QSW, as per :ref:`step 3 <demoral>` of the demoralisation procedure.

    :param G: Adjacency matrix, :math:`G`.
    :type G: :math:`N \\times N` SciPy CSR matrix

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :param perm_i: Use the :math:`i^\\text{th}` permutation of the Fourier matrices to orthogonalise the vertex subspaces of :math:`L^D`.
    :type perm_i: integer, optional

    :returns: A demoralised global Lindblad operator, :math:`L^D`.
    :rtype: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR matrix
    """

    row_ind = []
    col_ind = []
    data = []

    sub_counts = np.zeros((len(vsets),len(vsets)), dtype = np.complex128)
    for i in range(G.shape[0]):
        for j in range(G.shape[0]):
            if G[i,j] != 0:
                for k in range(len(vsets[i])):
                    for l in range(len(vsets[j])):
                        sub_counts[i,j] += G[i,j]

    for i in range(G.shape[0]):
        f = __fourier_matrix(len(vsets[i]), perm_i)
        for l in range(len(vsets[i])):
            index = 0
            for j in range(G.shape[0]):
                if G[i,j] != 0:
                    for k in range(len(vsets[j])):
                        row_ind.append(vsets[j][k])
                        col_ind.append(vsets[i][l])
                        data.append((1.0/np.sqrt(sub_counts[i,j])) * f[l, index])
                    index += 1

    shape = (vsets[-1][-1] + 1, vsets[-1][-1] + 1)

    return sp.sparse.csr_matrix((data,(col_ind, row_ind)), shape = shape)

def nm_H_loc(vsets):
    """
    Returns the locally rotating Hamiltonian, :math:`H_\\text{loc}`, as per :ref:`step 4 <demoral>` of the graph demoralised procedure.

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :returns: The locally rotating Hamiltonian, :math:`H_\\text{loc}`.
    :rtype: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR matrix
    """

    row_ind = []
    col_ind = []
    data = []

    for v in vsets:
        for i in range(len(v)):
            for j in range(len(v)):
                if v[i] != v[j]:
                    if (v[i] < v[j]) and (np.abs((v[i] - v[j])) == 1):
                        row_ind.append(v[i])
                        col_ind.append(v[j])
                        data.append(I)
                if (v[i] > v[j]) and (np.abs((v[i] - v[j])) == 1):
                        row_ind.append(v[i])
                        col_ind.append(v[j])
                        data.append(-I)

    shape = (vsets[-1][-1] + 1, vsets[-1][-1] + 1)

    return sp.sparse.csr_matrix((data,(row_ind, col_ind)), shape = shape)

def nm_rho_map(probs, vsets):
    """
    Maps the vertex populations, :math:`p(t_0)`, of a QSW on graph :math:`G`, to the corresponding vertex subspaces of the demoralised graph, :math:`G^D` (see Equation :eq:`eq:nm_rho_map`).

    :param probs: Vertex populations of :math:`\\rho(t_0)`.
    :type probs: :math:`1 \\times N` float array

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :returns: Density operator, :math:`\\rho^D(t_0)`
    :rtype: :math:`\\tilde{N} \\times \\tilde{N}` SciPy CSR matrix
    """

    rho = np.zeros(vsets[-1][-1] + 1)

    for i, prob in enumerate(probs):
        if prob > 0:
            rho[vsets[i]] = np.real(prob) / np.real(len(vsets[i]))

    return np.diag(rho)

def nm_measure(rho, vsets):
    """
    Maps the vertex populations of :math:`\\rho^D(t)` or :math:`(\\rho^D(t_0),...,\\rho^D(t_q))` to :math:`p(t)` or :math:`(p(t_0),...,p(t_q))` (see Equation :eq:`eq:nm_rho_map`).

    :param rho: :math:`\\rho^D(t)` or :math:`(\\rho^D(t_0),...,\\rho^D(t_q))`.
    :type rhp: :math:`\\tilde{N} \\times \\tilde{N}` or :math:`q \\times \\tilde{N} \\times \\tilde{N}` NumPy complex matrix

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :returns: :math:`p(t)` or :math:`(p(t_0),...,p(t_q))`.
    :rtype: :math:`N \\times N` or :math:`q \\times N \\times N` NumPy complex array.
    """

    if rho.ndim == 2:
        probs = []
        for v in vsets:
            probs.append(np.sum(rho[v,v]))
        return np.real(np.array(probs))
    else:
        populations = np.empty(
                (len(vsets),rho.shape[0]),
                dtype = np.float64)
        for i in range(rho.shape[0]):
            for j, v in enumerate(vsets):
                populations[j,i] = (np.sum(np.real(rho[i,:,:][v,v])))
        return populations

def nm_pop_inv_map(populations, vsets):
    """
    Maps the vertex populations of :math:`\\rho^D(t)` to :math:`\\rho(t)` (see Equation :eq:`eq:nm_rho_map`) given :math:`p^D(t)` or :math:`(p^D(t_0),...,p^D(t_q))`.

    :param rho: Vertex populations, :math:`p^D(t)` or :math:`(p^D(t_0),...,p^D(t_q))`.
    :type rhp: :math:`\\tilde{N} \\times \\tilde{N}` or :math:`q \\times \\tilde{N} \\times \\tilde{N}` NumPy complex array 

    :param vsets: Vertex subspaces, :math:`V^D`, generated via :meth:`~qsw_mpi.operators.nm_vsets`.
    :type vsets: integer arrays

    :returns: :math:`p(t)` or :math:`(p(t_0),...,p(t_q))`.
    :rtype: :math:`N \\times N` or :math:`q \\times N \\times N` NumPy complex array.
    """


    if populations.ndim == 1:

        nm_populations = np.empty(len(vsets), dtype = np.complex128)

        for i, v in enumerate(vsets):
            nm_populations[i] = np.sum(populations[v])

        return nm_populations

    else:

        nm_populations = np.empty(
                (len(vsets), populations.shape[1]),
                dtype = np.complex128)

        for j in range(populations.shape[1]):
            for i, v in enumerate(vsets):
                nm_populations[i,j] = np.sum(populations[:,j][v])

        return nm_populations
