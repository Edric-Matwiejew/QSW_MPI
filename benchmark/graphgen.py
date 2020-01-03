import numpy as np
import networkx as nx
import scipy.sparse as sparse
from pathlib import Path

print('Creating test graphs...')

Path("graphs/line_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/grid_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/random_graphs").mkdir(parents = True, exist_ok = True)
Path("graphs/complete_graphs").mkdir(parents = True, exist_ok = True)

np.random.seed(1)

r1 = 2
r2 = 21

r1_complete = 2
r2_complete =  10 #22

for i in range(r1,r2):

    graph = nx.path_graph(i**2, create_using=None)
    G = nx.to_scipy_sparse_matrix(graph)
    G = sparse.csr_matrix(G, dtype = np.float64)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.normal())

    sparse.save_npz('graphs/line_graphs/line_graph' + '_' + str(i**2) + '.npz', G)

    print('line graph' + str(i**2) + ' vertices')

for i in range(r1,r2):

    graph = nx.grid_graph(dim=[i,i])
    G = nx.to_scipy_sparse_matrix(graph)
    G = sparse.csr_matrix(G, dtype = np.float64)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.normal())


    sparse.save_npz('graphs/grid_graphs/grid_graph' + '_' + str(i) + '.npz', G)

    print('grid graph' + str(i**2) + ' vertices')

for i in range(r1,r2):

    graph = nx.gnm_random_graph(i**2, int((i**2)*np.log(i**2)), seed = None, directed = False)
    G = nx.to_scipy_sparse_matrix(graph)
    G = sparse.csr_matrix(G, dtype = np.float64)

    for j in range(0, len(G.data)):
            G.data[j] = G.data[j]*abs(np.random.normal())


    sparse.save_npz('graphs/random_graphs/random_graph' + '_' + str(i**2) + '.npz', G)

    print('random graph' + str(i**2) + ' vertices')

for i in range(r1_complete,r2_complete):

    graph = np.random.rand(i**2, i**2)
    np.fill_diagonal(graph,0)

    G = sparse.csr_matrix(graph)

    sparse.save_npz('graphs/complete_graphs/complete_graph' + '_' + str(i**2) + '.npz', G)

    print('complete graph' + str(i**2) + ' vertices')

print("Graph creation done.")
