#!/bin/bash

if [ -e QSW_MPI_Results ]; then
	rm -r QSW_MPI_Results
	rm -r plots
fi

mkdir -p QSW_MPI_Results
mkdir -p plots

export OMP_NUM_THREADS=1

echo "Started QSW_MPI series benchmark..."

mpiexec -N 4 python3 QSW_MPI_series.py $total_time

echo "done."

echo "Plotting results..."

python3 plot_line.py
python3 plot_complete.py

echo "done."
