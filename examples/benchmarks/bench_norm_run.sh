#!/bin/bash
export OMP_NUM_THREADS=1

mkdir graphs/line_graphs/sym
mkdir graphs/grid_graphs/sym
mkdir graphs/random_graphs/sym
mkdir graphs/complete_graphs/sym

echo "Line Graphs" >> benchmarks.dat
for file in graphs/line_graphs/*.npz
do
	time mpiexec -N 1 python benchmark_norms.py $file 'benchmark_lines_norms.csv'
done

echo "Square Grids" >> benchmarks.dat
for file in graphs/grid_graphs/*.npz
do
	time mpiexec -N 1 python benchmark_norms.py $file 'benchmark_grids_norms.csv'
done

echo "Random Graphs" >> benchmarks.dat
for file in graphs/random_graphs/*.npz
do
	time mpiexec -N 1 python benchmark_norms.py $file 'benchmark_randoms_norms.csv'
done

echo "Complete Graphs" >> benchmarks.dat
for file in graphs/complete_graphs/*.npz
do
	time mpiexec -N 1 python benchmark_norms.py $file 'benchmark_completes_norms.csv'
done





