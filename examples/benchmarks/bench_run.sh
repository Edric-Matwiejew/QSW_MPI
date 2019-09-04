#!/bin/bash
export OMP_NUM_THREADS=1

mkdir graphs/line_graphs/sym
mkdir graphs/grid_graphs/sym
mkdir graphs/random_graphs/sym
mkdir graphs/complete_graphs/sym

echo "FreeQSW benchmarks" > benchmarks.dat

echo "Line Graphs" >> benchmarks.dat
for file in graphs/line_graphs/*.npz
do
	(time mpiexec -N 4 python3 benchmark_steps.py $file 'benchmark_lines.csv') &>> benchmarks.dat
done

echo "Square Grids" >> benchmarks.dat
for file in graphs/grid_graphs/*.npz
do
	(time mpiexec -N 4 python3 benchmark_steps.py $file 'benchmark_grids.csv') &>> benchmarks.dat
done

echo "Random Graphs" >> benchmarks.dat
for file in graphs/random_graphs/*.npz
do
	(time mpiexec -N 4 python3 benchmark_steps.py $file 'benchmark_randoms.csv') &>> benchmarks.dat
done

echo "Complete Graphs" >> benchmarks.dat
for file in graphs/complete_graphs/*.npz
do
	(time mpiexec -N 4 python3 benchmark_steps.py $file 'benchmark_completes.csv') &>> benchmarks.dat
done





