import numpy as np
import networkx as nx
import scipy.sparse as sparse
from scipy.io import mmwrite
from pathlib import Path
from qsw_mpi.operators import symmetrise

print('Creating test graphs...')

Path("graphs/sym").mkdir(parents = True, exist_ok = True)

np.random.seed(1)

i = 5050
graph = nx.path_graph(i, create_using=None)
G = nx.to_scipy_sparse_matrix(graph)
G = sparse.csr_matrix(G, dtype = np.complex128)

for j in range(0, len(G.data)):
        G.data[j] = G.data[j]*abs(np.random.random())

G_sym = symmetrise(G)

mmwrite('graphs/line_graph' + '_' + str(i), G, field = 'complex', precision = 15)
mmwrite('graphs/sym/line_graph' + '_' + str(i) + '_sym', G_sym, field = 'complex', precision = 15)


i = 62
graph = nx.grid_graph(dim=[i,i])
G = nx.to_scipy_sparse_matrix(graph)
G = sparse.csr_matrix(G, dtype = np.complex128)

for j in range(0, len(G.data)):
        G.data[j] = G.data[j]*abs(np.random.random())

G_sym = symmetrise(G)

mmwrite('graphs/grid_graph' + '_' + str(i), G, field = 'complex', precision = 15)
mmwrite('graphs/sym/grid_graph' + '_' + str(i) + '_sym', G_sym, field = 'complex', precision = 15)

i = 2020
graph = nx.gnm_random_graph(i, int((i)*np.log(i)), seed = None, directed = False)
G = nx.to_scipy_sparse_matrix(graph)
G = sparse.csr_matrix(G, dtype = np.complex128)

for j in range(0, len(G.data)):
        G.data[j] = G.data[j]*abs(np.random.random())

G_sym = symmetrise(G)

mmwrite('graphs/random_graph' + '_' + str(i), G, field = 'complex', precision = 15)
mmwrite('graphs/sym/random_graph' + '_' + str(i) + '_sym', G_sym, field = 'complex', precision = 15)

i = 400
graph = np.random.rand(i, i)
np.fill_diagonal(graph,0)

G = sparse.csr_matrix(graph, dtype = np.complex128)

G_sym = symmetrise(G)

mmwrite('graphs/complete_graph' + '_' + str(i), G, field = 'complex', precision = 15)
mmwrite('graphs/sym/complete_graph' + '_' + str(i) + '_sym', G_sym, field = 'complex', precision = 15)

print("Graph creation done.")
