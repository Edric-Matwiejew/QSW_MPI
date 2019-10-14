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

mpi_comm = MPI.COMM_WORLD
rank = mpi_comm.rank

log = 'Application_2.csv'

files = ['graphs/A1.mtx', 'graphs/A2.mtx', 'graphs/A3.mtx', 'graphs/D1.mtx', 'graphs/D2.mtx', 'graphs/D3.mtx']

log = open(log, 'w')
log.write('graph, omega, integrated_populations, integrated_coherences, final_sink_pop, time, t2_population \n')

init_site = [1,3,7,2,6,14,3,9,21,4,12,28]

for j in range(len(files)):

    name = os.path.basename(files[j])[:-4]

    G = sp.csr_matrix(mmread(files[j]))

    gamma = 1
    H = qsw.operators.graph(gamma, G)
    L = qsw.operators.site_lindblads(H)
    print(H.todense())

    sinks = (np.array([0]), np.array([1]))

    omegas = np.array([0.01, 1])

    walk = qsw.MPI.walk(omegas[0], H, L, mpi_comm, sinks = sinks)
    walk.File('Application_2/' + name,"w")

    rho_0 = np.zeros(H.shape[0], dtype=np.complex128)
    rho_0[1:] = 1/float(H.shape[0])
    walk.initial_state(rho_0)

    t1 = 0
    t2 = 7000
    steps = 2**13

    for omega in omegas:

        walk.set_omega(omega)

        rho_t_series = walk.series(t1, t2, steps, target = 0, save = True, name = name + '_' + str(omega) + 'even')

        if rank == 0:

            if not os.path.isdir(r'plots/' + name):
                os.mkdir('plots/' + name)

            if omega == 0:
                graph = nx.from_scipy_sparse_matrix(G)
                qsw.plot.graph(graph, sinks = sinks, legend = False, sink_color='red', graph_color='orange')
                plt.savefig('plots/' + name + '/' + name + ".jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)

            h = (t2 - t1)/steps

            populations = qsw.measure.populations(rho = rho_t_series)

            transport_time = 0
            sink_pop_1 = 0
            it = enumerate(populations)
            next(it)
            for i, pop in enumerate(populations):
                sink_pop_2 = pop[H.shape[0]]
                if (((sink_pop_2 - sink_pop_1)/h) < (10**(-12))) and (sink_pop_1 != 0):
                    transport_time = i*h
                    break
                sink_pop_1 = sink_pop_2
            if transport_time == 0:
                print('Transport not completed over given time interval for ' + name)

            integrated_sites = np.empty(G.shape[0])
            for i in range(H.shape[0]):
                integrated_sites[i] = romb(populations[:,i], h)
            integrated_sites_sum = np.sum(integrated_sites)

            node_pairs, coherences = qsw.measure.coherences(rho = rho_t_series)

            integrated_coherences = np.empty(coherences.shape[1])
            for i in range(coherences.shape[1]):
                integrated_coherences[i] = romb(coherences[:,i], h)
            integrated_coherences_sum = np.sum(integrated_coherences)

            qsw.plot.population_lines(populations, [t1, t2], labels = True)
            plt.savefig('plots/' + name + '/' + 'pop' + '_' + str(omega) + '.jpeg')
            qsw.plot.coherence_lines(node_pairs, coherences, [t1, t2], labels = True)
            plt.savefig('plots/' + name + '/' + 'coh' + '_' + str(omega) + '.jpeg')

            log.write('{}, {}, {}, {}, {}, {}, {}\n'.format(
                files[j], omega, integrated_sites_sum, integrated_coherences_sum, populations[-1, H.shape[0]], \
                        transport_time, np.sum(populations[-1])))

            log.flush()


log.close()


