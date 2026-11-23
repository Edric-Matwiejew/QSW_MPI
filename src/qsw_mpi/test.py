import qsw_mpi.operators as op
import networkx as nx
import numpy as np
from scipy import sparse

graph = np.random.rand(3, 3)
np.fill_diagonal(graph,0)

G = sparse.csr_matrix(graph, dtype = np.complex128)

G_sym = op.sym(G)

vsets = op.nm_vsets(G)

H = op.nm_H(G,vsets)
L = op.nm_L(G,vsets)
H_loc = op.nm_H_loc(vsets)

print(vsets)
print(H_loc.todense())
