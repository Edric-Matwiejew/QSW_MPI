import qsw_mpi as qsw
import numpy as np
from glob import glob
from natsort import natsorted
from scipy.io import mmread
from scipy.sparse import csr_matrix as csr
import h5py
from mpi4py import MPI
from os.path import basename, exists
from time import time
import sys

def benchmark(omega, t1, t2, steps, files, files_sym, results, max_sim_time):

    total_sim_time = 0.0

    for f, fs in zip(files, files_sym):
        print(fs,flush=True)
        G = csr(mmread(fs))
        L = csr(mmread(f))

        sim_time = time()

        rho = np.zeros(G.shape)
        rho[0,0] = 1.0
        QSW = qsw.MPI.LQSW(omega, G, L, comm)
        QSW.initial_state(rho)
        QSW.series(t1, t2, steps)

        rho_t_series = QSW.gather_result()

        construct_time = comm.reduce(QSW.construct_time, MPI.MAX,0)
        rec_time = comm.reduce(QSW.rec_time, MPI.MAX, 0)
        one_norms_time = comm.reduce(QSW.one_norms_time, MPI.MAX, 0)
        series_time = comm.reduce(QSW.series_time, MPI.MAX, 0)
        so_nnz = comm.reduce(QSW.SO_col_indexes.shape[0], MPI.SUM,0)

        if rank == 0:
            step_time = 0
            diff_max_real = []
            diff_max_imag = []
            diff_min_real = []
            diff_min_imag = []
        for i in range(steps):
            QSW.step(t1 + i*(t2-t1)/steps)

            rho_t = QSW.gather_result()
            step_time_temp = comm.reduce(QSW.step_time, MPI.MAX, 0)
            if rank == 0:
                step_time += step_time_temp
                diff_max_real.append(np.max(np.abs(np.real(rho_t_series[i] - rho_t))))
                diff_max_imag.append(np.max(np.abs(np.imag(rho_t_series[i] - rho_t))))

        sim_time = time() - sim_time
        sim_time = comm.allreduce(sim_time, MPI.MAX)
        total_sim_time += sim_time

        if rank == 0:

            np.save('QSW_MPI_Results/diff_max_real_' + basename(f[:-4]),diff_max_real)
            np.save('QSW_MPI_Results/diff_max_imag_' + basename(f[:-4]),diff_max_imag)

            log = open(logname, 'a')
            log.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
            basename(f[:-4]),
            comm.Get_size(),
            G.shape[0],
            G.count_nonzero(),
            so_nnz,
            construct_time,
            rec_time,
            one_norms_time,
            series_time,
            step_time,
            sim_time))
            log.close()

        if total_sim_time > max_sim_time:
            break

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

omega = 0.1
t1 = 0.0
t2 = 100.0
steps = 1000
max_sim_time = np.float(sys.argv[1])

logname = 'QSW_MPI_Results/QSW_MPI_local_series_line.csv'
if rank == 0:
    if not exists(logname):
        log = open(logname, 'a')
        log.write('name,comm_size,dim,nnz,SO_nnz,SO_time,rec_time,one_norms_time,series_time,step_time,total_time\n')
        log.close()

files = natsorted(glob("../graphs/line_graphs/*.mtx"))
files_sym = natsorted(glob("../graphs/line_graphs/sym/*.mtx"))
benchmark(omega, t1, t2, steps, files, files_sym, "QSW_MPI_Results/line_graphs", max_sim_time)

logname = 'QSW_MPI_Results/QSW_MPI_local_series_complete.csv'
if rank == 0:
    if not exists(logname):
        log = open(logname, 'a')
        log.write('name,comm_size,dim,nnz,SO_nnz,SO_time,rec_time,one_norms_time,series_time,step_time,total_time\n')
        log.close()

files = natsorted(glob("../graphs/complete_graphs/*.mtx"))
files_sym = natsorted(glob("../graphs/complete_graphs/sym/*.mtx"))
benchmark(omega, t1, t2, steps, files, files_sym, "QSW_MPI_Results/complete_graphs", max_sim_time)
