import sys
import qsw_mpi
import numpy as np
from glob import glob
from scipy.io import mmread
from scipy.sparse import csr_matrix as csr
import h5py
from mpi4py import MPI
from os.path import basename, exists
from time import time
import resource

def benchmark(omega, t, files, files_sym, results, csv_name):

    total_sim_time = 0.0

    G = csr(mmread(files_sym))
    L = csr(mmread(files))

    sim_time = time()
    rho = np.identity(G.shape[0], dtype = np.complex128)/G.shape[0]
    QSW = qsw.MPI.LQSW(omega, G, L, comm)
    QSW.initial_state(rho)
    QSW.step(t)

    sim_time = time() - sim_time

    sim_time = comm.allreduce(sim_time, MPI.MAX)

    total_sim_time += sim_time

    rho_t = QSW.gather_result()

    construct_time = comm.reduce(QSW.construct_time, MPI.MAX,0)
    rec_time = comm.reduce(QSW.rec_time, MPI.MAX, 0)
    one_norms_time = comm.reduce(QSW.one_norms_time, MPI.MAX, 0)
    step_time = comm.reduce(QSW.step_time, MPI.MAX, 0)
    so_nnz = comm.reduce(QSW.SO_col_indexes.shape[0], MPI.SUM, 0)

    ave_comm_overhead = comm.reduce(np.sum(QSW.SO_num_send_inds), MPI.SUM, 0)

    peak_mem = comm.reduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, MPI.SUM,0)

    save_time = time()
    QSW.save_result(results, str(comm.Get_size()) + '/' + basename(files)[:-4], 'w')
    save_time = time() - save_time

    save_time = comm.reduce(save_time, MPI.MAX, 0)

    if rank == 0:

        peak_mem = peak_mem*(1024**(-2))

        log = open(csv_name, 'a')
        log.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
        basename(files[:-4]),
        comm.Get_size(),
        G.shape[0],
        G.count_nonzero(),
        so_nnz,
        construct_time,
        rec_time,
        one_norms_time,
        step_time,
        sim_time,
        peak_mem,
        save_time,
        ave_comm_overhead/float(comm.Get_size())))
        log.close()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

csv_name = sys.argv[4]

if rank == 0:
    if not exists(csv_name):
        log = open(csv_name, 'a')
        log.write('name,comm_size,dim,nnz,SO_nnz,SO_time,rec_time,one_norms_time,step_time,total_time,peak_memory,save_time,ave_comm_overhead\n')
        log.close()

omega = 0.1
t = 100.0

files = sys.argv[1]
files_sym = sys.argv[2]

benchmark(omega, t, files, files_sym, sys.argv[3], csv_name)
