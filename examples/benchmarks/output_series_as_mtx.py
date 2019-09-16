import sys
sys.path.append('../..')
import os
import numpy as np
from scipy import sparse as sp
from scipy.io import mmwrite
from scipy.sparse.linalg import expm_multiply
from mpi4py import MPI
import freeqsw as qsw
import h5py

def output_series_to_mtx(fi, mm_real, mm_imag, log):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        if os.path.exists(log):
            log = open(log, 'a')
        else:
            log = open(log, 'w')
            log.write('t, delta_re, delta_im\n')


    G = sp.csr_matrix(sp.load_npz(fi).todense(),dtype=np.complex128)
    H = qsw.operators.graph(1.0, G)
    L = qsw.operators.site_lindblads(H)
    qsw.operators.symmetrise(G)
    H = qsw.operators.graph(1.0, G)

    test_system = qsw.MPI.walk(0.1, H, L, comm)

    test_system.initial_state('even')

    rhot = test_system.series(0, 10, 1000, target = 0, precision = "sp")

    if rank == 0:

        f = h5py.File(mm_real)
        a_group_key = list(f.keys())[0]
        data = np.array(list(f[a_group_key]),dtype=np.float64)

        g = h5py.File(mm_imag)
        g_group_key = list(g.keys())[0]
        data_im = np.array(list(g[g_group_key]),dtype=np.float64)

        data = data + data_im*1j
        delta = np.array([np.flatten(data[i]-rhot[i]) for i in range(0,100)])
        delta_re = np.real(delta)
        delta_im = np.imag(delta)

        t = np.arange(0,100.01,0.01,dtype=np.float32)


mathematica_comparision(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
