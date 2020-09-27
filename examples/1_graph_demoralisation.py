"""
Here we provide a demonstration of graph moralisation
and the demoralisation correction scheme. This begins
by loading the required module; qsw_mpi, NumPy, CSR
matrix support from SciPy and MPI support from mpi4py.
"""
import qsw_mpi
import numpy as np
from scipy.sparse import csr_matrix as csr
from mpi4py import MPI
import matplotlib.pyplot as plt
import matplotlib

# Matplotlib parameters.
matplotlib.use("Agg")
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})
"""
As this system is small, only 1 MPI node will be used.
However, initialisation of an MPI communicator is
required to use the qsw_mpi.MPI module.
"""
comm = MPI.COMM_WORLD
"""
Below directed and undirected graphs are defined.
Here this is done by writing them directly into
the CSR format, where the arguments of 'csr_matrix' are
(non-zeros values,(row indices, column indices )) and graph
dimensions. The structure of the directed
graph and its undirected counterpart are as follows:

"""
if comm.Get_rank() == 0:
    """
    Note that code intended for 1 MPI process is nested inside
    a condinitional statement which checks for the target MPI rank.
    """

    print("G: \n 1 ---> 3 <--- 2")
    print("GU: \n 1 <--> 3 <--> 2")

G = csr(([1,1],([2,2],[0,1])), (3,3))
GU = csr(([1,1,1,1],([0,1,2,2],[2,2,0,1])), (3,3))

"""
First we will examine the behaviour of a G-QSW.
The Lindblad operator is defined as the directed adjacency matrix G and
The Hamilton as the undirected graph adjacency matrix GU.
Note that the Lindblad operator is stored in an array.
"""
Ls = [G]
H = GU
"""
Next we specify the starting state of the system as a pure state
at vertex 1. This may be achieved by either specifying the
density matrix completely, or by giving a list of probabilities,
in which case the off-diagonal entries are assumed to be 0. Here,
we have chosen the latter approach.
"""
rho_0 = np.array([1,0,0])
"""
The pieces are now in place for the definition of a walk object.
For a G-QSW the GQSW class is used.
"""
omega = 1.0
GQSW = qsw_mpi.MPI.GQSW(omega, H, Ls, comm)
"""
Finally, the initial state of the system is passed to the walk object.
"""
GQSW.initial_state(rho_0)
"""
Calculation of the system evolution at a single time point is done
though the 'step' method. After a simulation, we collect and matricise
the result at the root node using 'gather_result'.
"""
GQSW.step(t = 100)
rhot = GQSW.gather_result(root = 0)

if comm.Get_rank() == 0:
    """
    After the period of evolution we find that there is a non-zero
    probability of there being a walker at vertex 2, despite it
    having an in-degree of 0. This occurs due to there being a
    non-zero transition probability between vertices 1 and 3
    due to them having a common 'child' node, a phenomena termed
    'spontaneous moralisation'.
    """
    print("Global Interaction QSW, omega = 1: \n", np.real(rhot.diagonal()), flush = True)
"""
We will now demonstrate how to use QSW_MPI to apply the demoralisation
correction scheme. First we create the set of vertex subspaces.
"""
vsets = qsw_mpi.operators.nm_vsets(GU)
"""
These are then used with adjacency matrices G and GU to create
the Hamiltonian, Lindblad operators and rotating Hamiltonian
which capture the structure of the demoralised graph.
"""
H_nm = qsw_mpi.operators.nm_G(GU,vsets)
L_nm = [qsw_mpi.operators.nm_L(G,vsets)]
H_loc = qsw_mpi.operators.nm_H_loc(vsets)
"""
When creating the walk object, we called it with additional arguments
specifying the vertex subspaces and rotating Hamiltonian.
"""
nm_GQSW = qsw_mpi.MPI.GQSW(omega, H_nm, L_nm, comm, H_loc = H_loc, vsets = vsets)
"""
The initial system state is then mapped to the moralised graph as
and passed to the walk object. System propagation and
moralisation proceeds as before.
"""
rho_0_nm = qsw_mpi.operators.nm_rho_map(rho_0, vsets)
nm_GQSW.initial_state(rho_0)
nm_GQSW.step(t=100)
rhot_nm = nm_GQSW.gather_result()

if comm.Get_rank() == 0:
    """
    At t = 100 we find the system in a pure state at the sink node, as expected
    by the originating graph topology.

    >> [1.38389653e-87 0.00000000e+00 1.00000000e+00]
    """
    print("Non-moralisation GQSW, omega = 1: \n", qsw_mpi.operators.nm_measure(rhot_nm,vsets))

"""
As a further point of consideration we will now compare the dynamics of
the NM-G-QSW to an L-QSW on the same digraph, with Hamiltonian and
Lindblad operators defined as the adjacency matrices GU and G.
Note that the single-arc Lindblad operators are provided as a single CSR matrix.
"""

QSW = qsw_mpi.MPI.LQSW(omega, GU, G, comm)
QSW.initial_state(rho_0)

QSW.step(t = 100)

rho_t = QSW.gather_result()

if comm.Get_rank() == 0:
    """
    Evolving the state to $\rho(100)$ with \texttt{LQSW.step} yields,

    >> [-9.52705648e-18  0.00000000  1.00000000].

    Which corresponds to the state of the NM-G-QSW.
    """
    print("Local Interaction QSW, omega = 1: \n", np.real(rho_t.diagonal()), flush = True)

"""
The coherent evolution of the two systems is examined by first rebuilding 
the superoperators at omega = 0.
"""

nm_GQSW.set_omega(0)
QSW.set_omega(0)

nm_GQSW.step(t = 100)
QSW.step(t = 100)

nm_GQSW_rho_t = nm_GQSW.gather_result()
QSW_rho_t = QSW.gather_result()

if comm.Get_rank() == 0:
    """
    After which a step to $t = 100$ yields,

    >> [3.80773381e-07 9.98766244e-01 1.23337485e-03]

    for the NM-G-QSW and,

    >> [3.80773217e-07 9.98766244e-01 1.23337485e-03]

    for the L-QSW. In fact, for this particular system, the limiting dynamics
    of a NM-G-QSW correspond to that of a CTQW and CTRW, as is the case for the L-QSW.
    """
    print("Global Interaction Demoralised QSW, omega = 0: \n", qsw_mpi.operators.nm_measure(nm_GQSW_rho_t,vsets))
    print("Local Interaction QSW, omega = 0: \n", np.real(QSW_rho_t.diagonal()), flush = True)

"""
However, if we examine the time evolution of the two systems
at omega = 0.9 using the series method:
"""

nm_GQSW.set_omega(0.9)
QSW.set_omega(0.9)

nm_GQSW.series(t1=0,tq=25,steps=500)
QSW.series(t1=0,tq=25,steps=500)


anm_GQSW_rho_t = nm_GQSW.gather_result()
aQSW_rho_t = QSW.gather_populations()

if comm.Get_rank() == 0:
    """
    Notably different dynamics are observed. The nm-g-qsw results in a
    higher transfer of probability to the sink vertex and does not as
    readily decay to a quasi-stationary state.
    """
    plt.figure(figsize=(3.8,3.4))
    plt.xlabel('t', fontsize = 20)
    plt.ylabel(r'$p_{v_3}(t)$', fontsize = 20)
    plt.xticks(ticks=[0,250,500],labels=[0,12.5,25])
    plt.yticks(ticks=[0,0.5,1])
    plt.plot(qsw_mpi.operators.nm_measure(anm_GQSW_rho_t,vsets)[2,:])
    plt.plot(np.real(aQSW_rho_t[2,:]), '--')
    plt.savefig(
            '1_sink_dynamics',
            dpi=300,
            bbox_inches='tight',
            pad_inches = 0.05)
