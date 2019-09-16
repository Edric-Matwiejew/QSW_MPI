import sys
sys.path.append('../..')
import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import freeqsw as qsw
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

fi = 'graphs/line_graphs/line_graph_100.npz'

G = sp.csr_matrix(sp.load_npz(fi).todense(),dtype=np.complex128)
H = qsw.operators.graph(1.0, G)
L = qsw.operators.site_lindblads(H)
qsw.operators.symmetrise(G)
H = qsw.operators.graph(1.0, G)

test_system = qsw.MPI.walk(0.1, H, L, comm)

test_system.initial_state('even')

rhot = test_system.series(0, 10, 1000, target = 0, precision = "dp")

if rank == 0:

    f = h5py.File('../../../QSWalk/pt_line_re.hdf5')
    a_group_key = list(f.keys())[0]
    data = np.array(list(f[a_group_key]),dtype=np.float64)

    g = h5py.File('../../../QSWalk/pt_line_im.hdf5')
    g_group_key = list(g.keys())[0]
    data_im = np.array(list(g[g_group_key]),dtype=np.float64)

    data = data + data_im*1j
    delta = np.array([np.ndarray.flatten(data[i]-rhot[i]) for i in range(0,1000)])
    delta_re = np.real(delta)
    delta_im = np.imag(delta)
    t = np.arange(0,10,0.01,dtype=np.float32)

    print(np.max(np.abs(delta_re)))
    print(np.max(np.abs(delta_im)))

    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)
    plt.ylabel(r'$\Delta \rho$ $(1 \times 10^{-16})$')
    plt.xlabel('t')
    ax1.plot(t.T,delta_re, 'g+')
    ax1.plot(t.T, delta_im, 'b.')
    plt.tight_layout()
    plt.savefig("delta_line_dp.jpeg", bbox_inches = 'tight', pad_inches = 0.2)



