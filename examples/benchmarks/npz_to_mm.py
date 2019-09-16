import sys
sys.path.append('../..')
import os
import numpy as np
from scipy import sparse as sp
from scipy.io import mmwrite
import glob

files=glob.glob('./graphs/grid_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/grid_graphs/' + os.path.splitext(os.path.basename(fi))[-2], G)

files=glob.glob('./graphs/grid_graphs/sym/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/grid_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2], G)

#######################################################################################

files=glob.glob('./graphs/line_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/line_graphs/' + os.path.splitext(os.path.basename(fi))[-2], G)

files=glob.glob('./graphs/line_graphs/sym/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/line_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2], G)

########################################################################################

files=glob.glob('./graphs/random_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/random_graphs/' + os.path.splitext(os.path.basename(fi))[-2], G)

files=glob.glob('./graphs/random_graphs/sym/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/random_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2], G)

############################################################################################

files=glob.glob('./graphs/complete_graphs/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/complete_graphs/' + os.path.splitext(os.path.basename(fi))[-2], G)

files=glob.glob('./graphs/complete_graphs/sym/*npz')
sizes = []
for fi in files:
    sizes.append(glob.os.path.getsize(fi))
files = [x for _, x in sorted(zip(sizes,files))]

for fi in files:
    G = sp.load_npz(fi)
    mmwrite('graphs/complete_graphs/sym/' + os.path.splitext(os.path.basename(fi))[-2], G)












