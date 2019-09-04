import numpy as np
import networkx as nx
import scipy.sparse as sparse

for i in range(2,65):

        graph = nx.path_graph(i**2, create_using=None)
        G = nx.to_scipy_sparse_matrix(graph)
        G = sparse.csr_matrix(G, dtype = np.float64)

        for j in range(0, len(G.data)):
                G.data[j] = G.data[j]*abs(np.random.normal())

        sparse.save_npz('graphs/line_graphs/line_graph' + '_' + str(i**2) + '.npz', G)

        print('line ' + str(i))

        graph = nx.grid_graph(dim=[i,i])
        G = nx.to_scipy_sparse_matrix(graph)
        G = sparse.csr_matrix(G, dtype = np.float64)

        for j in range(0, len(G.data)):
                G.data[j] = G.data[j]*abs(np.random.normal())


        sparse.save_npz('graphs/grid_graphs/grid_graph' + '_' + str(i) + '.npz', G)

        print('grid ' + str(i))

        graph = nx.gnm_random_graph(i**2, int((i**2)*np.log(i**2)), seed = None, directed = False)
        G = nx.to_scipy_sparse_matrix(graph)
        G = sparse.csr_matrix(G, dtype = np.float64)

        for j in range(0, len(G.data)):
                G.data[j] = G.data[j]*abs(np.random.normal())


        sparse.save_npz('graphs/random_graphs/random_graph' + '_' + str(i**2) + '.npz', G)

        print('random ' + str(i))

        graph = np.random.rand(i**2, i**2)
        np.fill_diagonal(graph,0)

        G = sparse.csr_matrix(graph)

        sparse.save_npz('graphs/complete_graphs/complete_graph' + '_' + str(i**2) + '.npz', G)

        print('complete ' + str(i))

