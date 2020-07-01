import numpy as np
import networkx as nx
import scipy.sparse as sparse
from scipy.io import mmwrite
from pathlib import Path
from qsw_mpi.operators import symmetrise

print('Creating test graphs...')

Path("graphs/line_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/grid_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/random_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/complete_graphs").mkdir(parents = True, exist_ok = True)

Path("graphs/line_graphs/sym").mkdir(parents = True, exist_ok = True)
Path("graphs/grid_graphs/sym").mkdir(parents = True, exist_ok = True)
Path("graphs/random_graphs/sym").mkdir(parents = True, exist_ok = True)
Path("graphs/complete_graphs/sym").mkdir(parents = True, exist_ok = True)


np.random.seed(1)

r1 = 2
r2 = 56

r1_complete = 2
r2_complete =  21

for i in range(r1,r2):

    graph = nx.path_graph(i**2, create_using=None)
    G = nx.to_scipy_sparse_matrix(graph)
    G = sparse.csr_matrix(G, dtype = np.complex128)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.random())

    G_sym = symmetrise(G)

    mmwrite('graphs/line_graphs/line_graph' + '_' + str(i**2), G, field = 'complex', precision = 15)
    mmwrite('graphs/line_graphs/sym/line_graph' + '_' + str(i**2) + '_sym', G_sym, field = 'complex', precision = 15)

    print('line graph ' + str(i**2) + ' vertices')

for i in range(r1,r2):

    graph = nx.grid_graph(dim=[i,i])
    G = nx.to_scipy_sparse_matrix(graph)
    G = sparse.csr_matrix(G, dtype = np.complex128)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.random())

    G_sym = symmetrise(G)

    mmwrite('graphs/grid_graphs/grid_graph' + '_' + str(i), G, field = 'complex', precision = 15)
    mmwrite('graphs/grid_graphs/sym/grid_graph' + '_' + str(i) + '_sym', G_sym, field = 'complex', precision = 15)

    print('grid graph ' + str(i**2) + ' vertices')

for i in range(r1,r2):

    graph = nx.gnm_random_graph(i**2, int((i**2)*np.log(i**2)), seed = None, directed = False)
    rcm = list(nx.utils.rcm.cuthill_mckee_ordering(graph))
    G = nx.to_scipy_sparse_matrix(graph, nodelist = rcm)
    G = sparse.csr_matrix(G, dtype = np.complex128)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.random())

    G_sym = symmetrise(G)

    mmwrite('graphs/random_graphs/random_graph' + '_' + str(i**2), G, field = 'complex', precision = 15)
    mmwrite('graphs/random_graphs/sym/random_graph' + '_' + str(i**2) + '_sym', G_sym, field = 'complex', precision = 15)

    print('random graph ' + str(i**2) + ' vertices')

for i in range(r1_complete,r2_complete):

    graph = np.random.rand(i**2, i**2)
    np.fill_diagonal(graph,0)

    G = sparse.csr_matrix(graph, dtype = np.complex128)

    G_sym = symmetrise(G)

    mmwrite('graphs/complete_graphs/complete_graph' + '_' + str(i**2), G, field = 'complex', precision = 15)
    mmwrite('graphs/complete_graphs/sym/complete_graph' + '_' + str(i**2) + '_sym', G_sym, field = 'complex', precision = 15)

    print('complete graph ' + str(i**2) + ' vertices')

print("Graph creation done.")
