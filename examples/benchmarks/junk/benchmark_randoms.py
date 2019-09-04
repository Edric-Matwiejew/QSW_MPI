import sys
sys.path.append('../..')
import math
import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import networkx as nx
import freeqsw as qsw
import time
import glob
from memory_profiler import profile

@profile
def benchmark(files):
    for fi in files:
        total_time_start = time.time()

        if rank == 0:
            print(fi)

        G = sp.load_npz(fi)
        H = qsw.operators.graph(1.0, G)
        L = qsw.operators.site_lindblads(H)
        qsw.operators.symmetrise(H)

        test_system = qsw.MPI.walk(0.01, H, L, comm)
        test_system.initial_state('even')
        start = time.time()
        rhot = test_system.step(100, target = 0, precision = "sp")
        finish = time.time()

        total_time_end = time.time()

        if rank is 0:
            print(fi + 'step time: ' + str(finish - start))
            print(fi + 'total time: ' + str(total_time_end - total_time_start))

            pops = qsw.measure.populations(rho = rhot)

            print(fi + ' ' + str(np.sum(pops)))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

files=glob.glob('./graphs/random_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

benchmark(files)
