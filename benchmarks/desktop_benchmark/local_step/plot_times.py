import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

df_jl = pd.read_csv("Results/QSWalk_jl_local_step.csv")
df_m = pd.read_csv("Results/QSWalk_m_local_step.csv")
df_MPI = pd.read_csv("Results/QSW_MPI_local_step.csv")

df_MPI_1 = df_MPI[df_MPI["comm_size"] == 1]
max_comm = np.max(df_MPI["comm_size"].values)
df_MPI_max = df_MPI[df_MPI["comm_size"] == max_comm]

def plot_times(graph_type, savename):
    plt.figure(figsize=(5,5))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["total_time"], '.', label = "QSWalk.jl")
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["total_time"],'+', label = "QSWalk.m")
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["total_time"],'x', label = "QSW\_MPI (1)" )
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["total_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")")
    plt.xlabel(r'$\tilde{\mathcal{L}}$ Non-Zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.legend()
    plt.savefig(savename, dpi = 300)
    plt.close()

plot_times('line', 'Plots/line_total_times')
plot_times('grid', 'Plots/grid_total_times')
plot_times('random', 'Plots/random_total_times')
plot_times('complete', 'Plots/complete_total_times')

def exp_times(graph_type, savename):
    plt.figure(figsize=(5,5))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["step_time"], '.', label = "QSWalk.jl")
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["step_time"], '+', label = "QSWalk.m")
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["step_time"],'x', label = "QSW\_MPI (1)" )
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["step_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")")
    plt.xlabel(r'$\tilde{\mathcal{L}}$ Non-Zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.legend()
    plt.savefig(savename, dpi = 300)
    plt.close()

exp_times('line', 'Plots/line_exp_times')
exp_times('grid', 'Plots/grid_exp_times')
exp_times('random', 'Plots/random_exp_times')
exp_times('complete', 'Plots/complete_exp_times')

def construct_times(graph_type, savename):
    plt.figure(figsize=(5,5))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["SO_time"], '.', label = "QSWalk.jl")
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["SO_time"], '+', label = "QSWalk.m")
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_time"] \
            + df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["one_norms_time"],'x', label = "QSW\_MPI (1)" )
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_time"] \
            + df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["one_norms_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")")
    plt.xlabel(r'$\tilde{\mathcal{L}}$ Non-Zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.legend()
    plt.savefig(savename, dpi = 300)
    plt.close()

construct_times('line', 'Plots/line_construct_times')
construct_times('grid', 'Plots/grid_construct_times')
construct_times('random', 'Plots/random_construct_times')
construct_times('complete', 'Plots/complete_construct_times')
