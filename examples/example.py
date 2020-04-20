import sys
sys.path.append('../')

from mpi4py import MPI
import qsw_mpi as qsw
import numpy as np
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

matplot

mpi_comm = MPI.COMM_WORLD
rank = mpi_comm.Get_rank()

vertices = 4
Graph = nx.wheel_graph(vertices)
G = nx.to_scipy_sparse_matrix(Graph, dtype = np.complex128)

if rank == 0:
    matplotlib.use("Agg")
    ax = plt.subplot(label = 'graph')
    ax.axis('off')
    nx.draw_networkx(Graph, with_labels = True, node_color='orange')
    plt.savefig("wheel_graph.jpeg")

np.random.seed(1)
for i in range(G.count_nonzero()):
    G.data[i] = G.data[i]*np.random.random()

gamma = 1
H = qsw.operators.transition(gamma, G)
L = qsw.operators.lindblad(H)
qsw.operators.symmetrise(H)

omega = 0.001
wheel_graph = qsw.MPI.walk(omega, H, L, mpi_comm)

rho_0 = np.zeros(G.shape)
rho_0[0,0] = 1
wheel_graph.initial_state(rho_0)

t = 5
rho_t = wheel_graph.step(5, target = 0)

t1 = 0
t2 = 5

steps = 50

rho_t_series = wheel_graph.series(t1, t2, steps, target=0)

if rank == 0:

    pop_step = qsw.measure.populations(rho=rho_t)
    pop_series = qsw.measure.populations(rho=rho_t_series)

    print(pop_step)
    print(pop_series[50])
    print(np.sum(pop_step))

    vertex_pairs, cohs = qsw.measure.coherences(rho=rho_t_series)

    ax = plt.subplot(label = '2D plots')

    qsw.plot.population_lines(
            ax,
            pop_series,
            [t1, t2],
            labels = True)

    plt.savefig("population_lines.jpeg")

    ax.clear()

    qsw.plot.coherence_lines(
            ax,
            vertex_pairs,
            cohs,
            [t1, t2],
            labels = True)

    plt.savefig("coherence_lines.jpeg")

    ax.clear()

wheel_graph.initial_state('even')

wheel_graph.set_omega(0.5)

sources = ([0], [0.7])
sinks = ([3], [0.8])

if rank == 0:

    qsw.plot.graph(
            Graph,
            sources = sources,
            sinks = sinks)

    plt.savefig("augmented_wheel_graph.jpeg")

wheel_graph_augmented = qsw.MPI.walk(
        omega,
        H,
        L,
        mpi_comm,
        sources = sources,
        sinks = sinks)

wheel_graph_augmented.File('usage_example', action = "w")

wheel_graph_augmented.initial_state('sources')

wheel_graph_augmented.series(
        t1,
        t2,
        steps,
        save = True,
        name = "series 2",
        chunk_size = 5)

if rank == 0:

    file = qsw.io.File('usage_example')

    pop_series = qsw.measure.populations(File = file, series_name = "series 2")

    ax_3D = plt.subplot(projection = '3d')

    qsw.plot.population_bars( ax_3D, pop_series, [t1, t2])

    plt.savefig("population_bars.jpeg")

    ax_3D.clear()

    vertex_pairs, cohs = qsw.measure.coherences(File = file, series_name = "series 2")

    qsw.plot.coherence_bars(
            ax_3D,
            vertex_pairs,
            cohs,
            [t1, t2])

    plt.savefig("coherence_bars.jpeg")


wheel_graph_augmented = qsw.MPI.load_walk(omega, 'usage_example', mpi_comm)

if rank == 0:
    file.list_series()
