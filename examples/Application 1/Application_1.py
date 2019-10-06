import sys
sys.path.append('../../')
from mpi4py import MPI
import freeqsw as qsw
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import mmread
import glob
import os
from scipy import sparse as sp
from scipy.integrate import romb
import networkx as nx
from scipy.io import mmwrite

mpi_comm = MPI.COMM_WORLD
rank = mpi_comm.rank

log = 'Application_1.csv'

files=glob.glob('isomorphic_graphs/*.mtx')

log = open(log, 'w')
log.write('graph, omega, integrated_populations, integrated_coherences, final_sink_pop, time, t2_population \n')

for j in range(len(files)):

    name = os.path.basename(files[j])[:-4]

    G = sp.csr_matrix(mmread(files[j]))

    gamma = 1
    H = qsw.operators.graph(gamma, G)
    L = qsw.operators.site_lindblads(H)
    print(H.todense())

    sinks = (np.array([2]), np.array([3]))

#    omegas = np.arange(0, 1.01, 0.01)
#    print(omegas)
    omegas = [0.0]

    walk = qsw.MPI.walk(omegas[0], H, L, mpi_comm, sinks = sinks)
    print(walk.M_values)
    print(walk.M_col_indexes)
    print(walk.M_row_starts)
    mmwrite(name, sp.csr_matrix((walk.M_values, walk.M_col_indexes - 1, walk.M_row_starts - 1), shape = (64,64)))
    walk.File('results/' + name,"w")

    rho_0 = np.zeros(H.shape, dtype=np.complex128)
    rho_0[3,3] = 1
    walk.initial_state(rho_0)

    t1 = 0
    t2 = 50000
    steps = 2**16

#    for omega in omegas:
#
#        walk.set_omega(omega)
#
#        rho_t_series = walk.series(t1, t2, steps, target = 0, save = True, name = name + '_' + str(omega))
#
#        if rank == 0:
#
#            if not os.path.isdir(r'plots/' + name):
#                os.mkdir('plots/' + name)
#
#            if omega == 0:
#                graph = nx.from_scipy_sparse_matrix(G)
#                qsw.plot.graph(graph, sinks = sinks, legend = False, sink_color='red', graph_color='orange')
#                plt.savefig('plots/' + name + '/' + name + ".jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)
#
#            h = (t2 - t1)/steps
#
#            populations = qsw.measure.populations(rho = rho_t_series)
#
#            transport_time = 0
#            sink_pop_1 = 0
#            it = enumerate(populations)
#            next(it)
#            for i, pop in enumerate(populations):
#                sink_pop_2 = pop[7]
#                if (((sink_pop_2 - sink_pop_1)/h) < (10**(-12))) and (sink_pop_1 != 0):
#                    transport_time = i*h
#                    break
#                sink_pop_1 = sink_pop_2
#            if transport_time == 0:
#                print('Transport not completed over given time interval for ' + name)
#
#            integrated_sites = np.empty(G.shape[0])
#            for i in range(7):
#                integrated_sites[i] = romb(populations[:,i], h)
#            integrated_sites_sum = np.sum(integrated_sites)
#
#            node_pairs, coherences = qsw.measure.coherences(rho = rho_t_series)
#
#            integrated_coherences = np.empty(coherences.shape[1])
#            for i in range(coherences.shape[1]):
#                integrated_coherences[i] = romb(coherences[:,i], h)
#            integrated_coherences_sum = np.sum(integrated_coherences)
#
#            #qsw.plot.population_lines(populations, t1, t2, labels = True)
#           # plt.savefig('plots/' + name + '/' + 'pop' + '_' + str(omega) + '.jpeg')
#            #qsw.plot.coherence_lines(node_pairs, coherences, t1, t2, labels = True)
#            #plt.savefig('plots/' + name + '/' + 'coh' + '_' + str(omega) + '.jpeg')
#
#            log.write('{}, {}, {}, {}, {}, {}, {}\n'.format(
#                files[j], omega, integrated_sites_sum, integrated_coherences_sum, populations[-1,7], \
#                        transport_time, np.sum(populations[-1])))
#
#            log.flush()
#
#
#log.close()
#
#
