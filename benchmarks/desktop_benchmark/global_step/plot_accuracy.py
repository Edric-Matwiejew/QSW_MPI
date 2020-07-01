import h5py
import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 20})
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

        m_mean_cmplx = []
        m_max_cmplx = []
        m_min_cmplx = []
        m_mean_real = []
        m_max_real = []
        m_min_real = []

        for m, mpi in zip(m_keys, comm_keys):
            try:
                rho_m = np.array(df_m[m]).view(np.complex128)
                rho_mpi = np.array(df_mpi[str(i) + '/' + mpi]).view(np.complex128)

                m_mean_cmplx.append(np.mean(np.imag(rho_mpi - rho_m)))
                m_max_cmplx.append(np.max(np.imag(rho_mpi - rho_m)))
                m_min_cmplx.append(np.min(np.imag(rho_mpi - rho_m)))
                m_mean_real.append(np.mean(np.real(rho_mpi - rho_m)))
                m_max_real.append(np.max(np.real(rho_mpi - rho_m)))
                m_min_real.append(np.min(np.real(rho_mpi - rho_m)))
            except:
                break

        validation_data[str(i)] = {
            'm_mean_cmplx':m_mean_cmplx,
            'm_max_cmplx':m_max_cmplx,
            'm_min_cmplx':m_min_cmplx,
            'm_mean_real':m_mean_real,
            'm_max_real':m_max_real,
            'm_min_real':m_min_real
        }

    return validation_data

def dif_plot(valdat, prefix, color, comm_size, comp, plot_name):
    xs = [i**2 for i in range(1,len(valdat[comm_size][prefix + '_mean_' + comp]) + 1)]
    plt.figure(figsize = (5,5))
    plt.plot(xs, valdat[comm_size][prefix + '_mean_' + comp], '1', color = color)
    plt.plot(xs, valdat[comm_size][prefix + '_max_' + comp], '.', color = color)
    plt.plot(xs, valdat[comm_size][prefix + '_min_' + comp], '.', color = color)
    plt.xlabel("Verticies")
    plt.ylabel(r"min$(\Delta\rho(t))$    \hspace{1.8cm}   max($\Delta\rho(t)$)")
    plt.yscale('symlog', linthreshy = 1e-16)
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ylim(-1e-9,1e-9)
    plt.hlines(0, 1, max(xs) + 1,linestyles= '--', color = 'grey')
    plt.savefig(plot_name, dpi = 300)
    plt.close()

linedat = compare(df_m_line, df_mpi_line, 'line')
griddat = compare(df_m_grid, df_mpi_grid, 'grid')
randomdat = compare(df_m_random, df_mpi_random, 'random')
completedat = compare(df_m_complete, df_mpi_complete, 'complete')

for key in linedat.keys():
    dif_plot(linedat, 'm', 'green', key, 'real', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_real_dif.jpg')
    dif_plot(linedat, 'm', 'blue', key, 'cmplx', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_imag_dif.jpg')

    dif_plot(griddat, 'm', 'green', key, 'real', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_grid_real_dif.jpg')
    dif_plot(griddat, 'm', 'blue', key, 'cmplx', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_grid_imag_dif.jpg')

    dif_plot(randomdat, 'm', 'green', key, 'real', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_random_real_dif.jpg')
    dif_plot(randomdat, 'm', 'blue', key, 'cmplx', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_random_imag_dif.jpg')

    dif_plot(completedat, 'm', 'green', key, 'real', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_complete_dif.jpg')
    dif_plot(completedat, 'm', 'blue', key, 'cmplx', 'Plots/' + key + '_processes_' + 'QSW_MPI_vs_QSWalk_m_line_complete_dif.jpg')



