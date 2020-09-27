import h5py
import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

df_m_line = h5py.File("Results/line_graphs_QSWalk_m.h5", 'r')
df_mpi_line = h5py.File("Results/line_graphs_QSW_MPI.h5", 'r')

df_m_grid = h5py.File("Results/grid_graphs_QSWalk_m.h5", 'r')
df_mpi_grid = h5py.File("Results/grid_graphs_QSW_MPI.h5", 'r')

df_m_random = h5py.File("Results/random_graphs_QSWalk_m.h5", 'r')
df_mpi_random = h5py.File("Results/random_graphs_QSW_MPI.h5", 'r')

df_m_complete = h5py.File("Results/complete_graphs_QSWalk_m.h5", 'r')
df_mpi_complete = h5py.File("Results/complete_graphs_QSW_MPI.h5", 'r')

def compare(df_m, df_mpi, graph_name):

    m_keys = natsorted(list(df_m.keys()))

    mpi_keys = []
    for key in df_mpi.keys():
        mpi_keys.append(natsorted(list(df_mpi[key].keys())))

    validation_data = {}
    for i, comm_keys in zip(np.array(list(df_mpi.keys()),dtype = int), mpi_keys):

        m_max_cmplx = []
        m_mean_real = []
        m_max_real = []

        for m, mpi in zip(m_keys, comm_keys):
            try:
                rho_mpi = np.array(df_mpi[str(i) + '/' + mpi]).view(np.complex128)
                rho_m = np.array(df_m[m]).view(np.complex128)

                m_max_cmplx.append(np.abs(np.max(np.imag(rho_mpi - rho_m))))
                m_max_real.append(np.abs(np.max(np.real(rho_mpi - rho_m))))
            except:
                break

        validation_data[str(i)] = {
            'm_max_cmplx':m_max_cmplx,
            'm_max_real':m_max_real,
        }

    return validation_data

def dif_plot(valdat, prefix, comm_size, plot_name):
    xs = [i**2 for i in range(1,len(valdat[comm_size][prefix + '_max_real']) + 1)]
    plt.figure(figsize = (3.8,2.8))
    plt.plot(xs, valdat[comm_size][prefix + '_max_real'], 'x', color = 'green', markersize=8, label = 'Re')
    plt.plot(xs, valdat[comm_size][prefix + '_max_cmplx'], '+', color = 'blue', markersize = 10, label = 'Im')
    plt.xlabel("verticies")
    plt.ylabel(r"max$(|\Delta\rho(t)|)$")
    plt.yscale('log')
    plt.ylim(None,10e-9)
    x_max = np.max(xs)
    plt.xlim(None,x_max)
    plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
    plt.yticks([10e-16,10e-14,10e-12,10e-10])
    plt.hlines(0, 1, max(xs) + 1,linestyles= '--', color = 'grey')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.legend()
    plt.savefig(plot_name, dpi = 300, bbox_inches = 'tight', edgecolour='none', pad_inches = 0.05)
    plt.close()

linedat = compare(df_m_line, df_mpi_line, 'line')
griddat = compare(df_m_grid, df_mpi_grid, 'grid')
randomdat = compare(df_m_random, df_mpi_random, 'random')
completedat = compare(df_m_complete, df_mpi_complete, 'complete')

for key in linedat.keys():
    dif_plot(linedat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_dif.jpg')
    dif_plot(linedat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_dif.jpg')

    dif_plot(griddat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_grid_dif.jpg')
    dif_plot(griddat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_grid_dif.jpg')

    dif_plot(randomdat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_random_dif.jpg')
    dif_plot(randomdat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_random_dif.jpg')

    dif_plot(completedat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_complete_dif.jpg')
    dif_plot(completedat, 'm', key, 'Plots/local_' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_complete_dif.jpg')
