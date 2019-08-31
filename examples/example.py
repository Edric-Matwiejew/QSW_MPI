import sys
sys.path.append('../')
import math
import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import networkx as nx
import freeqsw as qsw

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#Graph = nx.star_graph(10)
Graph = nx.path_graph(10)
G = nx.to_scipy_sparse_matrix(Graph, dtype=np.complex128)

np.random.seed(10)
for i in range(G.count_nonzero()):
    G.data[i] = G.data[i]*np.random.normal()

H = qsw.operators.graph(1.0, G)
L = qsw.operators.site_lindblads(H)
qsw.operators.symmetrise(H)

source_sites = np.array([30,70])
source_rates = np.array([1.0,2.0])
sink_sites = np.array([0, 99])
sink_rates = np.array([1, 0.5])

sources = (source_sites, source_rates)
sinks = (sink_sites, sink_rates)

test_system = qsw.MPI.walk(0.00, H, L, comm, sources = sources, sinks = sinks)
#test_system = qsw.MPI.walk(1.00, H, L, comm)

test_system.file('test', action = "w")
#test_system.initial_state('sources')

#test_system = qsw.MPI.load_walk(0.1, 'test', comm)
test_system.initial_state('sources')
test_system.step(100, save = True, precision = "sp")
test_system.series(0.0, 30.0, 500, save = True)
test_system.set_omega(0.1)
test_system.series(0.0, 30.0, 500, save = True)

if rank is 0:

    walks = qsw.io.File('test.qsw')

    pops_1 = qsw.measure.populations(File = walks, series_name = 'series 1')
    pops_2 = qsw.measure.populations(File = walks, series_name = 'series 2')

    print(pops_1[3])
    print(np.sum(pops_1[3]))
    print(pops_2[3])
    print(np.sum(pops_2[3]))

    qsw.plot.animate(Graph, pops_1, 0, 100.0, 'MyAnimation', sources = sources, sinks = sinks, node_size = 1000, framerate = 15, use_labels = False)
    #qsw.plot.animate(Graph, pops, 0, 100.0, 'MyAnimation', node_size = 4000, framerate = 15)

    walks.close()
