import sys
sys.path.append('../..')
import os
import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import freeqsw as qsw
import time
from memory_profiler import memory_usage

def benchmark(fi, log):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        if os.path.exists(log):
            log = open(log, 'a')
        else:
            log = open(log, 'w')
            log.write('file, system size, M nnz, norm_diff, construct t, reconcile t, one norms t, series t, total t, memory\n')


    total_start = time.time()

    G = sp.load_npz(fi)
    H = qsw.operators.graph(1.0, G)
    L = qsw.operators.site_lindblads(H)

    qsw.operators.symmetrise(H)

    test_system = qsw.MPI.walk(0.1, H, L, comm)

    test_system.initial_state('even')

    step_start = time.time()
    rhot = test_system.series(0, 10, 1000, target = 0, precision = "dp")
    step_end = time.time()

    total_end = time.time()

    if rank == 0:
        pops = qsw.measure.populations(rho = rhot)

    memory = memory_usage()[0]

    total_memory = comm.reduce(memory, MPI.SUM)

    if rank == 0:
        norm_diff = np.max(np.abs(1.0 - np.sum(pops,axis=1)))
        system_size = H.shape[0]
        M_nnz = test_system.M_values.shape[0]
        step_time = step_end - step_start
        total_time = total_end - total_start

        log.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
            fi, system_size, M_nnz, norm_diff, test_system.construction_time, \
                    test_system.reconcile_time, test_system.one_norm_time, \
                        step_time, total_time, total_memory))

        log.close()

benchmark(sys.argv[1], sys.argv[2])
