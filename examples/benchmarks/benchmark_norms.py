import sys
sys.path.append('../..')
import os
import numpy as np
from scipy import sparse as sp
from mpi4py import MPI
import freeqsw as qsw
import time
import glob
from memory_profiler import memory_usage

def benchmark(fi, log):
    
    files=glob.glob(fi +'*.npz')
    sizes = []
    for file in files:
        sizes.append(glob.os.path.getsize(file))
    files = [x for _, x in sorted(zip(sizes,files))]

    print(files)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    with open(log, 'w') as f:
    
        for file in files:

            G = sp.load_npz(file)
            H = qsw.operators.graph(1.0, G)
            L = qsw.operators.site_lindblads(H)
          
            qsw.operators.symmetrise(H)

            test_system = qsw.MPI.walk(0.1, H, L, comm)

            M = sp.csr_matrix((test_system.M_values, test_system.M_col_indexes - 1, test_system.M_row_starts - 1),(test_system.M_rows, test_system.M_rows))

            norm_diff = []
            for i in range(0, test_system.p):
                norm = np.max(np.sum(np.abs(M.todense()**(i+1)),1))
                norm_diff.append(100*(norm - test_system.one_norms[i])/norm)
            print(norm_diff)

            f.write(str(file) + ',' + str(norm_diff) + '\n')
            f.flush()
            
benchmark(sys.argv[1], sys.argv[2])
