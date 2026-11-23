"""Dicke superradiance from 32 initially excited atoms in a common reservoir."""

import numpy as np
from mpi4py import MPI
from scipy.sparse import csr_matrix, diags

from qsw_mpi.MPI import GKSL

comm = MPI.COMM_WORLD
atoms = 32
# Basis state k is the symmetric Dicke state with k excited atoms.
k = np.arange(atoms + 1)
rates = k * (atoms - k + 1)
jump = diags(np.sqrt(rates[1:]), offsets=1, format="csr", dtype=complex)
hamiltonian = csr_matrix(jump.shape, dtype=complex)
open_system = GKSL(hamiltonian, [jump], comm)

initial = np.zeros(atoms + 1)
initial[-1] = 1.0
open_system.initial_state(initial)
times = np.linspace(0.0, 0.5, 101)
open_system.series(times[0], times[-1], len(times) - 1)
populations = open_system.gather_populations()

if comm.rank == 0:
    # Each collective jump emits one photon.
    emission = rates @ populations.real
    peak = np.argmax(emission)
    print("Time    Photon emission rate")
    for sample in (0, peak, len(times) - 1):
        print(f"{times[sample]:.3f}   {emission[sample]:.3f}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 16})
    fig, ax = plt.subplots(figsize=(3.8, 3.4), layout="tight")
    ax.plot(times, emission)
    ax.set_xlabel(r"Time $t$", fontsize=20)
    ax.set_ylabel("Photon emission rate", fontsize=20)
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(bottom=0)
    fig.savefig("gksl_superradiance.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
