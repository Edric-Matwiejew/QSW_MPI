import numpy as np
from scipy import sparse as sp

def populations(rho = None, File = None, step_name = None, series_name = None):
    """Return the populations of a density matrix or density matrix series.
    """

    if rho is not None:
        if rho.ndim is 2:
            site_pops = np.real(np.diag(rho))
            return site_pops
        elif rho.ndim is 3:
            site_pops = np.array([np.real(np.diag(rho_step)) for rho_step in rho])
            return site_pops

    if File is not None:
        if step_name is not None:
            site_pops = np.real(np.diag(np.array(File.File['steps'][str(step_name)], dtype = np.complex128)))
            return site_pops

        elif series_name is not None:
            steps = File.File['series'][str(series_name)].attrs['steps'] + 1
            site_pops = np.array([np.real(np.diag(np.array(File.File['series'][str(series_name)][i,:,:], \
                    dtype = np.complex128))) for i in range(steps)])
            return site_pops

def coherences(rho = None, File = None, step_name = None, series_name = None):
    """Return the inter-site coherences of a density matrix or density matrix series.
    """

    if rho is not None:
        if rho.ndim is 2:

            node_pairs = np.triu_indices_from(rho, k=1)
            cohs = np.asarray(abs(rho[node_pairs]))

            return node_pairs, cohs

        elif rho.ndim is 3:

            node_pairs = np.triu_indices_from(rho[0], k=1)
            cohs = np.asarray(abs(rho[0][node_pairs]))
            it = iter(rho)
            next(it,None)

            for step in it:
                cohs = np.vstack((cohs,np.asarray(abs(step[node_pairs]))))

            return node_pairs, cohs

    if File is not None:
        if step_name is not None:

            rho = np.array(File.File['steps'][str(step_name)], dtype = np.complex128)

            node_pairs = np.triu_indices_from(rho[0], k=1)
            cohs = np.asarray(abs(rho[0][node_pairs]))

            return node_pairs, cohs

        elif series_name is not None:

            steps = File.File['series'][str(series_name)].attrs['steps'] + 1

            rho = np.array(File.File['series'][str(series_name)][0], dtype = np.complex128)

            node_pairs = np.triu_indices_from(rho[0], k=1)
            cohs = np.asarray(abs(rho[0][node_pairs]))

            it = iter(range(steps))
            next(it, None)

            for i in it:
                cohs = np.vstack((cohs,np.asarray(abs(np.array(File.File['series'][str(series_name)][i,:,:], dtype = np.complex128)[node_paris]))))

            return node_pairs, cohs
