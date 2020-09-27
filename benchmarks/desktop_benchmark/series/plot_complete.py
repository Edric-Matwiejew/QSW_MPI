import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from natsort import natsorted
import colour
import pandas as pd
from matplotlib import use

use("Agg")

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'figure.autolayout': True})

diffs_max_real = natsorted(glob('QSW_MPI_Results/*max_real_complete*.npy'))
diffs_max_imag = natsorted(glob('QSW_MPI_Results/*max_imag_complete*.npy'))

light_orange = colour.Color("#ffff99")
colors = list(light_orange.range_to(colour.Color("#db6a00"),len(diffs_max_real)))

plt.figure(figsize=(3.8,2.8))
for i, (re_diff, im_diff) in enumerate(zip(diffs_max_real,diffs_max_imag)):
    reals = np.load(re_diff, allow_pickle=True)
    imags = np.load(im_diff, allow_pickle=True)

    plt.plot(reals,'.', color = colors[i].rgb, alpha = 0.5, markeredgewidth = 0.0, label = 'Re')
    plt.plot(imags,'.', color = colors[i].rgb, alpha = 0.5, markeredgewidth = 0.0, label = 'Im')


plt.hlines(0, 1, 1000 + 1,linestyles= '--', color = 'grey')
plt.ylabel(r"max$(|\Delta \rho(t)|)$")
plt.xlabel("time step")
plt.yscale('log')
plt.ylim(1e-17,1e-13)
plt.yticks([10e-17,10e-16,10e-15,10e-14,10e-13])
plt.legend()
plt.savefig("plots/series_delta_complete.jpg", dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
plt.clf()


data_df = pd.read_csv('QSW_MPI_Results/QSW_MPI_local_series_complete.csv')

plt.figure(figsize=(3.8,3.4))
plt.plot(data_df['SO_nnz'], data_df['series_time'],'+',label="series",markersize=10)
plt.plot(data_df['SO_nnz'],data_df['step_time'],'x',label="step",markersize=8)
plt.plot(data_df['SO_nnz'],data_df['step_time'].values/1000,'1',label="step mean",markersize=10)
x_max = plt.xlim()[1]
y_max = plt.ylim()[1]
plt.xticks(np.array([0,round((x_max/2)/10**np.floor(np.log10(x_max/2)),1)*10**np.floor(np.log10(x_max/2)), round(x_max/10**np.floor(np.log10(x_max)),1)*10**np.floor(np.log10(x_max))]))
plt.yticks(np.array([0,round((y_max/2)/10**np.floor(np.log10(y_max/2)),1)*10**np.floor(np.log10(y_max/2)), round(y_max/10**np.floor(np.log10(y_max)),1)*10**np.floor(np.log10(y_max))]))
plt.ylabel("time (s)")
plt.ticklabel_format(axis='x', style='sci', scilimits = (0,0))
plt.ticklabel_format(axis='y', style='sci', scilimits = (0,0))
plt.xlabel(r"$\tilde{\mathcal{L}}$ Non-zeros")
plt.legend()
plt.savefig("plots/step_vs_series_times_complete.jpg", dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)
plt.close()
