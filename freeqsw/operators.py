import numpy as np
from scipy import sparse as sp
import freeqsw.foperators as foperators

def check_indices(A):
    if not A.has_sorted_indices:
        A.sort_indices()

def graph(gamma, A):

    check_indices(A)

    G_indptr, G_indices, G_data, G_nnz = foperators.graph(gamma, A.indptr, A.indices, A.data)

    G = sp.csr_matrix((G_data[0:G_nnz], G_indices[0:G_nnz], G_indptr), A.shape)

    return G

def site_lindblads(G):

    check_indices(G)

    L_indptr, L_indices, L_data, L_nnz = foperators.site_lindblads(G.indptr, G.indices, G.data)

    L = sp.csr_matrix((L_data[0:L_nnz], L_indices[0:L_nnz], L_indptr), G.shape)

    return L

def symmetrise(G):

    check_indices(G)

    foperators.symmetrise(G.indptr, G.indices, G.data)

