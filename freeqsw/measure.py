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
