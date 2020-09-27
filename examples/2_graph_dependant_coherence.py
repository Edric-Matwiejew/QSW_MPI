import qsw_mpi
from scipy.sparse import csr_matrix
from mpi4py import MPI
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib

# Matplotlib parameters.
matplotlib.use("Agg")

comm = MPI.COMM_WORLD

"""
Steady state of an L-QSW on a 2-branching tree graph of depth 5:
"""
Graph = nx.balanced_tree(2,5)
G = nx.to_scipy_sparse_matrix(Graph)

H = qsw_mpi.operators.transition(1.0, G)

"""
Local-interction Lindblad operators are derived from the
cannonical markov chain transition matrix.
"""
M = qsw_mpi.operators.markov_chain(G)
L = qsw_mpi.operators.local_lindblads(M)

omega = 0.5
QSW = qsw_mpi.MPI.LQSW(omega, H, L, comm)

"""
The system begins in a maximally mixed state.
"""
QSW.initial_state('mixed')

QSW.step(t=100)
rho_t = QSW.gather_result()

# Plot the graph structure and steady state.
if comm.Get_rank() == 0:
    plt.figure(figsize=(4,4))
    plt.axis('off')
    plt.imshow(np.log(np.abs(rho_t)), cmap = 'Purples')
    plt.savefig('2_tree_state',dpi=300, bbox_inches='tight', pad_inches = 0.05)
    plt.close()

    pos = graphviz_layout(Graph, prog="twopi")
    plt.figure(figsize=(4,4))
    plt.axis('off')
    nx.draw(Graph,pos=pos,node_color='purple')
    plt.tight_layout()
    plt.savefig(
            '2_tree_graph',
            dpi=300,
            bbox_inches='tight',
            pad_inches = 0.05)
    plt.close()

"""
Steady state on an L-QSW on a cycle graph of 60 verticies.
"""
N = 63
Graph = nx.cycle_graph(N)

G = nx.to_scipy_sparse_matrix(Graph)
H = qsw_mpi.operators.transition(1.0, G)

M = qsw_mpi.operators.markov_chain(G)
L = qsw_mpi.operators.local_lindblads(M)

QSW = qsw_mpi.MPI.LQSW(omega, H, L, comm)

QSW.initial_state('mixed')

QSW.step(t=100)
rho_t = QSW.gather_result()

# Plot the graph structure and steady state.
if comm.Get_rank() == 0:
    plt.figure(figsize=(4,4))
    plt.axis('off')
    pos = graphviz_layout(Graph)
    nx.draw(Graph,pos=pos,node_color='green',node_size=100)
    plt.tight_layout()
    plt.savefig('2_cycle_graph',dpi=300, bbox_inches='tight', pad_inches = 0.05)
    plt.close()

    plt.figure(figsize=(4,4))
    plt.axis('off')
    colmap = matplotlib.cm.Greens
    colmap.set_under('white')

    rhot_plot = np.abs(np.log(np.abs(rho_t),
                out=np.zeros_like(np.abs(rho_t)),
                where=(np.abs(rho_t)!=0)))

    plt.imshow(rhot_plot, cmap = colmap,vmin=0.0001)
    plt.savefig(
            '2_cycle_state',
            dpi=300,
            bbox_inches='tight',
            pad_inches = 0.05)
    plt.close()
