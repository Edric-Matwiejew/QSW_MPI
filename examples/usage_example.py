# This example will simulate a QSW on a directed wheel graph, shown is Figure _.
# We first must import mpi4py and FreeQSW, these are the minimum requirements to run a FreeQSW simulation.
# Also used in this example are the modules networkx, numpy and matplotlib.
from mpi4py import MPI
import qsw_mpi as qsw
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# The next step is to step up the MPI enviroment. This involves creating an MPI communicator object and defining the process rank.

mpi_comm = MPI.COMM_WORLD
rank = mpi_comm.Get_rank()

# EXAMPLE 1

# A graph is then constructed at each node using the nx.wheel_graph(4), where 4 is the number of nodes. The networkx graph object is labeled with sequential integers to ensure that the numbering of its nodes remains consistent durring calls to the networkx library.  The graph then convereted to a Scipy sparse CSR matrix.

nodes = 4
Graph = nx.wheel_graph(nodes)
G = nx.to_scipy_sparse_matrix(Graph, dtype=np.complex128)

fig, ax = plt.subplots(1, figsize=(2.8, 2.2))
ax.axis('off')
nx.draw_networkx(Graph, with_labels = True, node_color='orange')
plt.tight_layout()
plt.savefig("wheel_graph.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)

# To create a directed graph we will assign random weights between [0, 1) to the edges of G. Note the use of np.random.seed(1), this ensures that each MPI process generates the same sequence of random numbers.

np.random.seed(1)
for i in range(G.count_nonzero()):
    G.data[i] = G.data[i]*np.random.random()

#From this  a directed graph laplacian is created with $\gamma = 1$. This is used to obtain the Lindblad matrix and then symmetrised as per Equation _  to obtain H.

gamma = 1
H = qsw.operators.graph(gamma, G)
L = qsw.operators.site_lindblads(H)
qsw.operators.symmetrise(H)

# We then instantiate a walk object, on creation it constructs the distributed superoperator $\LL$ and determines its 1-norm series. Below $\omega = 0.1$.

omega = 0.1
wheel_graph = qsw.MPI.walk(omega, H, L, mpi_comm)

# To perform a simulation we must define the initial state of $\rho(0)$. This is passed to the walk object which vectorises and partitions $\rho(0)$ over the MPI communicator. Here note that  $rho(0) = |1><1|$ is being specified by the indices (0,0).

rho_0 = np.zeros(G.shape, dtype=np.complex128)
rho_0[0,0] = 1
wheel_graph.initial_state(rho_0)

# A single time point is properageted through use of the 'step' method. Let's do this for $t=5$ and return the result to $rank = 0$.
t = 5
rho_t = wheel_graph.step(5, target = 0)

# We can obtain a time series using the 'series' method. Let's do that for t = [0, 5] with 50 steps. Again the result is being returned to $rank = 0$.
t1 = 0
t2 = 5
steps = 50
rho_t_series = wheel_graph.series(t1, t2, steps, target=0)

# To analyse the results we must take case that we only act of the arrays rho_t and rho_t_series from the target rank, as they will not be defined elsewhere. To do so, we simply check that the process $rank$ is equal to the $target$ rank. Analysis may proceed by acting directly on the rho_t and rho_t_series numpy arrays, or through use of the qsw.measurements sub-module. For example we can obtain the populations from rho_t and rho_t_series by:

if rank == 0:

    pop_step = qsw.measure.populations(rho=rho_t)
    pop_series = qsw.measure.populations(rho=rho_t_series)

# Where the last element of pop_series is equal to pops:

    print(pop_step)
    print(pop_series[50])

# And as expected the populations total to 1.

    print(np.sum(pop_step))

# Internode coherences can likewise be extracted with the coherences method.

    node_pairs, cohs = qsw.measure.coherences(rho=rho_t_series)

# These results can be visualised using the FreeQSW plot sub-module.

    #qsw.plot.population_lines(pop_series, [t1, t2], labels = True)
    #plt.savefig("population_lines.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)
    #qsw.plot.coherence_lines(node_pairs, cohs, [t1, t2], labels = True)
    #plt.savefig("coherence_lines.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)

# Which are shown in Figures _ to _.

# The 'wheel_graph' walk object we may carry out simulations with different inital conditions by redfining the inital state. A sample number of default inital states are support. For example we can define the initial state as being an equal superposition across all nodes by:

wheel_graph.initial_state('even')

# We can also perform walks using  `wheel_graph' for other values of $\omega$ by redefining this parameter:

wheel_graph.set_omega(0.5)

# EXAMPLE 2

# Let's now expand the simulation to include aborption and emission by attaching a `source' to site 0 and a `sink' to site `3' of Graph. To do this we first define sources and sink tuples. The each case the first element contains a list nodes of connecting  points, and the second a list of rates.

sources = (np.array([0]), np.array([0.7]))
sinks = (np.array([3]), np.array([0.8]))

# The sesulting graph can be visualized using the plot submodule:

#qsw.plot.graph(Graph, sources = sources, sinks = sinks, size = (4.5, 3.5), graph_color='orange', source_color='green', sink_color='red')
#plt.savefig("augmented_wheel_graph.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)

# By Equation _,  these additions result in a structrual change to H and L, it is therefore nessecary to intianstite  a new walk object.

wheel_graph_augmented = qsw.MPI.walk(   omega, H, L, mpi_comm, \
                                        sources = sources, \
                                        sinks = sinks)

# When running remote, or large, simulations it more convienent to save simulation results directly to disc.  To do so we first create a .qsw file to contain the results of walks carried out using the wheel_graph_augmented object and the data required to recpnstruct the system at a future date.

wheel_graph_augmented.File('usage_example', action = "w")

# We then create an initial state whereby the walker is distributed accorss the source node(s) at t = t1. Below a walk is again carried out for the interval $t = [0, 5]$ with 50 steps with the addition of a 'name' parameter with which uniquely identifies this walk in the useage_example.qsw. The chunk_size parameter controls the number of time-steps between each write to disk, this is useful for the case where gathering the entire series to one root process exceeds the avaliable memory.

wheel_graph_augmented.initial_state('sources')
rhot = wheel_graph_augmented.series(t1, t2, steps, save = True, \
                                    name = "series 1", \
                                    chunk_size = 5, target = 0)

# For analysis the node populations may be extracted directly from the compressed .qsw file, thus avoiding the need to load the entire series to memory. As a systems increase in size it is useful to make use of a another plotting dimension for greater clairity.

if rank == 0:
    file = qsw.io.File('usage_example')
    pop_series = qsw.measure.populations(rho=rhot)
    node_pairs, cohs = qsw.measure.coherences(rho=rho_t_series)

    #qsw.plot.population_bars(   pop_series, [t1, t2], \
    #                            t_tick_freq = 10, \
    #                            t_round = 1)
    #plt.savefig("population_bars.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)
    #plt.figure(figsize=(1,1))
    #qsw.plot.coherence_bars(node_pairs, cohs, [t1, t2], \
    #                        t_tick_freq = 10, t_round = 1)
    #plt.savefig("coherence_bars.jpeg", bbox_inches = 'tight', pad_inches = 0.2, dpi=300)

# At a later time, if we wish to conduct further simulations with the same system, it may be initialized using 'useage_example.qsw', which contains all of the arrays nessecary to reconstruct the same super-operator.

wheel_graph_augmented = qsw.MPI.load_walk(omega, 'usage_example', mpi_comm)
