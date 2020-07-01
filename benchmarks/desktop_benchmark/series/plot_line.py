import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from natsort import natsorted
import colour
import pandas as pd
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

diffs_max_real = natsorted(glob('QSW_MPI_Results/*max_real_line*.npy'))
diffs_max_imag = natsorted(glob('QSW_MPI_Results/*max_imag_line*.npy'))
diffs_min_real = natsorted(glob('QSW_MPI_Results/*min_real_line*.npy'))
diffs_min_imag = natsorted(glob('QSW_MPI_Results/*min_imag_line*.npy'))

green = colour.Color("#98fb98")
colors = list(green.range_to(colour.Color("#013220"),len(diffs_max_real)))
for i, diff in enumerate(diffs_min_real):
    plt.plot(np.load(diff, allow_pickle=True),'.', color = colors[i].rgb)
    plt.yscale('symlog', linthreshy = 1e-16)

for i, diff in enumerate(diffs_max_real):
    plt.plot(np.load(diff, allow_pickle=True),'.', color = colors[i].rgb)
    plt.yscale('symlog', linthreshy = 1e-16)

plt.hlines(0, 1, 1000 + 1,linestyles= '--', color = 'grey')
plt.ylabel(r"min($\Delta \rho(t)$) \hspace{1.5cm} max($\Delta \rho(t)$)")
plt.xlabel("time step")
plt.ylim(-1e-13,1e-13)
plt.savefig("plots/delta_real_line.jpg", dpi = 300)
plt.clf()

blue = colour.Color("#afeeee")
colors = list(blue.range_to(colour.Color("#00008b"),len(diffs_max_real)))
for i, diff in enumerate(diffs_min_imag):
    plt.plot(np.load(diff, allow_pickle=True),'.', color = colors[i].rgb)
    plt.yscale('symlog', linthreshy = 1e-16)

for i, diff in enumerate(diffs_max_imag):
    plt.plot(np.load(diff, allow_pickle=True),'.', color = colors[i].rgb)
    plt.yscale('symlog', linthreshy = 1e-16)

plt.hlines(0, 1, 1000 + 1,linestyles= '--', color = 'grey')
plt.ylabel(r"min($\Delta (\rho(t)$) \hspace{1.5cm} max($\Delta \rho(t)$)")
plt.xlabel("time step")
plt.ylim(-1e-13,1e-13)
plt.savefig("plots/delta_imag_line.jpg", dpi = 300)
plt.clf()

data_df = pd.read_csv('QSW_MPI_Results/QSW_MPI_local_series_line.csv')

plt.figure(figsize=(5,5))
plt.plot(data_df['SO_nnz'], data_df['series_time'],'+',label="series")
plt.plot(data_df['SO_nnz'],data_df['step_time'],'x',label="step")
plt.plot(data_df['SO_nnz'],data_df['step_time'].values/1000,'1',label="step mean")
plt.ylabel("time (s)")
plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
plt.xlabel(r"$\tilde{\mathcal{L}}$ Non-zeros")
plt.legend()
plt.savefig("plots/step_vs_series_times_line.jpg", dpi = 300)
plt.close()
