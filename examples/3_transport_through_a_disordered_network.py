"""
This example makes use of QSW_MPI's efficient time series calculation
to reproduce results which illustrate that transport through a
disordered network as a function of omega can be closely modelled as transport through
an energetically disordered dimer.
"""
import numpy as np
from scipy.sparse import csr_matrix
from mpi4py import MPI
import qsw_mpi
"""
For this we will use the plotting library Matplotlib, and least-squares
optimization and numerical integration methods provided by SciPy.
"""
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import least_squares as leastsq
from scipy import integrate

# Matplotlib parameters.
matplotlib.use("TkAgg")
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.autolayout': True})

"""
Specifically, we will be studying a system of N points randomly
distributed in a sphere of radius 1. To which a source will
be attached at point 1 and a sink at point N.
"""
def random_3D_points(N):
    points = np.empty((3,N), dtype = np.float64)
    for i in range(N):
        while True:
            points[:,i] = np.random.uniform(-1, 1, (3,))
            if np.linalg.norm(points[:,i]) < 1:
                break
    return points
"""
Which are undergoing dipole-dipole interactions:
"""
def dipole_interaction(points):

    H = np.empty((points.shape[1],points.shape[1]))

    for i in range(points.shape[1]):
        for j in range(points.shape[1]):
            distance = np.linalg.norm(points[:,i] - points[:,j])
            if distance > 0:
                H[i,j] = 1/distance**3
            else:
                H[i,j] = 0

    return H
"""
To study the efficiency of transport through the network we consider
the Expected Survival Time (EST):
"""
def expected_survival_time(omegas, H, L, sources, sinks, t1, t2, steps):

    if rank == 0:
        EST = []

    for omega in omegas:

        """
        An L-QSW is initialised the optional arguments 'sources'
        and 'sinks', which are tuples for the form ([vertex],[rate]).
        """
        walk = qsw_mpi.MPI.LQSW(omega, H, L, mpi_comm, sources = sources, sinks = sinks)
        """
        To create an initial state distributed equally across all attached
        sources we pass the string 'sources' to the 'initial_state' method.
        """
        walk.initial_state('sources')
        """
        We can then calculated the time series of the system evolution using
        the 'walk' method:
        """
        walk.series(t1,t2,steps)
        """
        And retrieve the vertex populations.
        """
        series_pops = walk.gather_populations()
        """
        With which the EST as a function of omega is calculated.
        """
        if rank == 0:
            exit_index = []
            for i, pop in enumerate(series_pops[:][-1]):
                exit_index.append(np.real(1.0 - pop))
            EST.append(integrate.simps(exit_index, dx = 0.1))

    if rank == 0:
        return np.array(EST)
    else:
        return None

"""
We can now calculate EST for the dipole system for omega = 0.1,0.2,...,1.
"""

np.random.seed(1)

mpi_comm = MPI.COMM_WORLD
rank = mpi_comm.Get_rank()

N = 10

points = random_3D_points(N)

H = csr_matrix(dipole_interaction(points))
L = qsw_mpi.operators.local_lindblads(H)

sources = ([0],[0.5])
sinks = ([N - 1],[0.5])

t1 = 0
t2 = 500
steps = 1000
omegas = np.linspace(0.1,1,10)

expected_survival_times = expected_survival_time(omegas, H, L, sources, sinks, t1, t2, steps)
"""
And broadcast the result from the root MPI nodes to all other in its
communicator.
"""
EST = mpi_comm.bcast(expected_survival_times, root = 0)
"""
To fit the response of EST as function of omega for the dimer to that of the
dipole network we will make use of least squares optimisation with the
following objective function:
"""
def objective(dimer_params):
    """
    Where we will optimise over parameters V, Gamma and gamma,
    """
    global network, dimer

    V = dimer_params[0]
    Gamma = dimer_params[1]
    gamma = dimer_params[2]

    H = csr_matrix([[0,-V],[-V, Delta]])
    L = qsw_mpi.operators.local_lindblads(H)
    sources = ([0],[Gamma])
    sinks = ([1],[gamma])

    expected_survival_time_dimer = expected_survival_time(omegas, H, L, sources, sinks, t1, t2, steps)

    EST_dimer = mpi_comm.bcast(expected_survival_time_dimer, root = 0)

    if rank == 0:
        """
        and progress of the optimisation is visualised using Matplotlib.
        """
        plt.clf()
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$\eta(\omega)$')
        plt.xlim(0,1.1)
        plt.xticks(ticks = [0,0.2,0.4,0.6,0.8,1.0])
        network = plt.plot(omegas, EST, 'o', label = "Network", markersize = 8)
        dimer = plt.plot(omegas, EST_dimer, '+', label = "Dimer", color = 'r', markersize = 10)
        plt.pause(0.01)

    return EST_dimer - EST

Delta = 1.5

V = 0.5
gamma = 0.5
Gamma = 0.5

plt.figure(figsize=(5,4))

result = leastsq(objective,[V,gamma,Gamma], bounds = (10e-13, np.inf), xtol = 10e-8, ftol = 10e-8)

if rank == 0:

    print(result)

    plt.legend()
    plt.savefig('3_dimer_fit', dpi = 300, bbox_inches='tight', pad_inches = 0.05)
    plt.close()
