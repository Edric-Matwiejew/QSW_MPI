import numpy as np
from scipy import sparse as sp
import qsw_mpi.foperators as foperators

def check_indices(A):
#CSR column index array must be sorted.
    if not A.has_sorted_indices:
        A.sort_indices()

def graph(gamma, A):
    """Generate the graph Laplacian from an adjacency matrix.

    :param gamma: Scaling factor, must satisty 0 <:math:`\gamma`.
    :type gamma: numpy.real64

    :param A: Graph adjacency matrix.
    :type A: csr_matrix

    :rtype: csr_matrix
    """
    if gamma < 0:
        raise ValueError("0 < gamma")

    check_indices(A)

    G_indptr, G_indices, G_data, G_nnz = foperators.graph(gamma, A.indptr, A.indices, A.data)

    G = sp.csr_matrix((G_data[0:G_nnz], G_indices[0:G_nnz], G_indptr), A.shape)

    return G

def site_lindblads(G):
    """Generate a Lindblad operator matrix from a graph Laplacician.

    :param A: Graph Laplacian.
    :type A: csr_matrix

    :rtype: csr_matrix
    """
    check_indices(G)

    L_indptr, L_indices, L_data, L_nnz = foperators.site_lindblads(G.indptr, G.indices, G.data)

    L = sp.csr_matrix((L_data[0:L_nnz], L_indices[0:L_nnz], L_indptr), G.shape)

    return L

def symmetrise(G):
    """Symmetrise a csr_matrix with a symmetric structure.

    G is altered so as to satisfy :math:`G_{ij} = max(G_{ij}, G_{ji})`.

    :param G: Matrix to symmetrise.
    :type G: csr_matrix

    .. warning::
        The matrix must be structurally symmetric in terms of its non-zero entires.
    """
    check_indices(G)

    foperators.symmetrise(G.indptr, G.indices, G.data)
