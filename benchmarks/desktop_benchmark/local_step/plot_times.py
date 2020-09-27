import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

df_jl = pd.read_csv("Results/QSWalk_jl_local_step.csv")
df_m = pd.read_csv("Results/QSWalk_m_local_step.csv")
df_MPI = pd.read_csv("Results/QSW_MPI_local_step.csv")

df_MPI_1 = df_MPI[df_MPI["comm_size"] == 1]
max_comm = np.max(df_MPI["comm_size"].values)
df_MPI_max = df_MPI[df_MPI["comm_size"] == max_comm]

def plot_times(graph_type, savename):
    plt.figure(figsize=(3.8,3.4))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["total_time"], '.', label = "QSWalk.jl",markersize=10)
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["total_time"],'+', label = "QSWalk.m",markersize=10)
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["total_time"],'x', label = "QSW\_MPI (1)" , markersize=8)
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["total_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")",markersize=10)
    x_max = plt.xlim()[1]
    y_max = plt.ylim()[1]
    plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
    plt.yticks(np.array([0,round((y_max/2)/10**np.floor(np.log10(y_max/2)),1)*10**np.floor(np.log10(y_max/2)), round(y_max/10**np.floor(np.log10(y_max)),1)*10**np.floor(np.log10(y_max))]))
    plt.xlabel(r'$\tilde{\mathcal{L}}$ non-zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.savefig(savename, dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
    plt.clf()

plot_times('line', 'Plots/local_line_total_times')
plot_times('grid', 'Plots/local_grid_total_times')
plot_times('random', 'Plots/local_random_total_times')
plot_times('complete', 'Plots/local_complete_total_times')

def exp_times(graph_type, savename):
    plt.figure(figsize=(3.8,3.4))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["step_time"], '.', label = "QSWalk.jl", markersize=10)
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["step_time"], '+', label = "QSWalk.m",markersize=10)
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["step_time"],'x', label = "QSW\_MPI (1)" ,markersize=8)
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["step_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")",markersize=10)
    plt.xlabel(r'$\tilde{\mathcal{L}}$ non-zeros')
    x_max = plt.xlim()[1]
    y_max = plt.ylim()[1]
    plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
    plt.yticks(np.array([0,round((y_max/2)/10**np.floor(np.log10(y_max/2)),1)*10**np.floor(np.log10(y_max/2)), round(y_max/10**np.floor(np.log10(y_max)),1)*10**np.floor(np.log10(y_max))]))
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.savefig(savename, dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
    plt.close()

exp_times('line', 'Plots/local_line_exp_times')
exp_times('grid', 'Plots/local_grid_exp_times')
exp_times('random', 'Plots/local_random_exp_times')
exp_times('complete', 'Plots/local_complete_exp_times')

def construct_times(graph_type, savename):
    plt.figure(figsize=(3.8,3.4))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["SO_time"], '.', label = "QSWalk.jl",markersize=10)
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["SO_time"], '+', label = "QSWalk.m",markersize=10)
    plt.plot(df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_time"],'x', label = "QSW\_MPI (1)" ,markersize=8)
    plt.plot(df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"], df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")",markersize=10)
    x_max = plt.xlim()[1]
    y_max = plt.ylim()[1]
    plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
    plt.yticks(np.array([0,round((y_max/2)/10**np.floor(np.log10(y_max/2)),1)*10**np.floor(np.log10(y_max/2)), round(y_max/10**np.floor(np.log10(y_max)),1)*10**np.floor(np.log10(y_max))]))
    plt.xlabel(r'$\tilde{\mathcal{L}}$ non-zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.savefig(savename, dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
    plt.close()

construct_times('line', 'Plots/local_line_construct_times')
construct_times('grid', 'Plots/local_grid_construct_times')
construct_times('random', 'Plots/local_random_construct_times')
construct_times('complete', 'Plots/local_complete_construct_times')

def exp_plus_norm_times(graph_type, savename):
    plt.figure(figsize=(3.8,3.4))
    plt.plot(df_jl[df_jl["name"].str.startswith(graph_type)]["SO_nnz"], df_jl[df_jl["name"].str.startswith(graph_type)]["step_time"], '.', label = "QSWalk.jl",markersize=10)
    plt.plot(df_m[df_m["name"].str.startswith(graph_type)]["SO_nnz"], df_m[df_m["name"].str.startswith(graph_type)]["step_time"], '+', label = "QSWalk.m",markersize=10)
    plt.plot(
        df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["SO_nnz"],
        df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["step_time"] \
        + df_MPI_1[df_MPI_1["name"].str.startswith(graph_type)]["one_norms_time"],'x', label = "QSW\_MPI (1)" ,markersize=8)
    plt.plot(
        df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["SO_nnz"],
        df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["step_time"] \
        + df_MPI_max[df_MPI_max["name"].str.startswith(graph_type)]["one_norms_time"],'1', label = "QSW\_MPI (" + str(max_comm) + ")",markersize=10)
    x_max = plt.xlim()[1]
    y_max = plt.ylim()[1]
    plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
    plt.yticks(np.array([0,round((y_max/2)/10**np.floor(np.log10(y_max/2)),1)*10**np.floor(np.log10(y_max/2)), round(y_max/10**np.floor(np.log10(y_max)),1)*10**np.floor(np.log10(y_max))]))
    plt.xlabel(r'$\tilde{\mathcal{L}}$ non-zeros')
    plt.ylabel('time (s)')
    plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
    plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
    plt.savefig(savename, dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
    plt.close()

exp_plus_norm_times('line', 'Plots/local_line_exp_plus_norm_times.jpg')
exp_plus_norm_times('grid', 'Plots/local_grid_exp_plus_norm_times.jpg')
exp_plus_norm_times('random', 'Plots/local_random_exp_plus_norm_times.jpg')
exp_plus_norm_times('complete', 'Plots/local_complete_exp_plus_norm_times.jpg')
