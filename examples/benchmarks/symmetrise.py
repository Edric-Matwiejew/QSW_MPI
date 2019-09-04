import sys
sys.path.append('../..')
import os
import numpy as np
from scipy import sparse as sp
import freeqsw as qsw
import glob

files=glob.glob('./graphs/grid_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    A = sp.csr_matrix(G, dtype=np.complex128)
    qsw.operators.symmetrise(A)
    sp.save_npz('graphs/grid_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2] + '_sym.npz', A)

files=glob.glob('./graphs/line_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    A = sp.csr_matrix(G, dtype=np.complex128)
    qsw.operators.symmetrise(A)
    sp.save_npz('graphs/line_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2] + '_sym.npz', A)

files=glob.glob('./graphs/random_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    A = sp.csr_matrix(G, dtype=np.complex128)
    qsw.operators.symmetrise(A)
    sp.save_npz('graphs/random_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2] + '_sym.npz', A)

files=glob.glob('./graphs/complete_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    A = sp.csr_matrix(G, dtype=np.complex128)
    qsw.operators.symmetrise(A)
    sp.save_npz('graphs/complete_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2] + '_sym.npz', A)





