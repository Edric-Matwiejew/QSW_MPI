import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import re

def natural_sort(l):

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

test_graphs = ["line", "grid", "random", "complete"]

plt.figure(figsize = (25, 5))

for i, graph in enumerate(test_graphs):
    results = natural_sort(glob("results/*" + graph + "*.csv"))

    dfs = []
    for result in results:
        dfs.append(pd.read_csv(result))

    plt.subplot(1, 4, i + 1)
    for j, df in enumerate(dfs):
        plt.scatter(df["Super-operator nnz"], df["total t"], label = results[j])

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.xlabel("Super-operator (non-zeros)")
    plt.ylabel("time (s)")

plt.savefig("Step Benchmark.jpeg")
