import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from natsort import natsorted
import colour
import pandas as pd
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

diffs_max_real = natsorted(glob('QSW_MPI_Results/*max_real_line*.npy'))
diffs_max_imag = natsorted(glob('QSW_MPI_Results/*max_imag_line*.npy'))

light_orange = colour.Color("#ffff99")
colors = list(light_orange.range_to(colour.Color("#db6a00"),len(diffs_max_real)))

plt.figure(figsize=(3.8,2.8))
for i, (re_diff, im_diff) in enumerate(zip(diffs_max_real,diffs_max_imag)):
    reals = np.load(re_diff, allow_pickle=True)
    imags = np.load(im_diff, allow_pickle=True)

    max_diffs = []
    for re, im in zip(reals, imags):
        max_diffs.append(np.max([re,im]))

    plt.plot(max_diffs,'.', color = colors[i].rgb)

plt.ylabel(r"max$(\mathrm{Re}|\Delta \rho(t)|,\mathrm{Im}|\Delta \rho(t)|)$")
plt.xlabel(r"time steps")
plt.yscale('log')
plt.ylim(1e-17,1e-13)
plt.yticks([10e-17,10e-16,10e-15,10e-14,10e-13])
plt.savefig("plots/series_delta_line.jpg", dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
plt.clf()

data_df = pd.read_csv('QSW_MPI_Results/QSW_MPI_local_series_line.csv')

plt.figure(figsize=(3.8,3.8))
plt.plot(data_df['SO_nnz'], data_df['series_time'],'+',label="series",markersize=10)
plt.plot(data_df['SO_nnz'],data_df['step_time'],'x',label="step",markersize=8)
plt.plot(data_df['SO_nnz'],data_df['step_time'].values/1000,'1',label="step mean",markersize=10)
plt.ylabel(r"time (s)")
plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
plt.xlabel(r"$\tilde{\mathcal{L}}$ Non-zeros")
plt.legend()
plt.savefig("plots/step_vs_series_times_line.jpg", dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
plt.close()
