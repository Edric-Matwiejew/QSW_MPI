#!/bin/bash

mkdir -p results

python3 graphgen.py

nodes=(1 2 3 4)

echo "Starting QSW_MPI benchmarks:"

for node in ${nodes[@]}; do

	terminal_out="results/benchmarks_$node.dat"

	echo "$node MPI nodes" > $terminal_out

	echo "Line Graphs" >> $terminal_out
	for file in `ls -v graphs/line_graphs/*.npz`; do
		echo "Running $file on $node nodes..."
		(time mpiexec -N $node python3 steps.py $file "results/benchmark_lines_$node.csv") &>> $terminal_out
	done

	echo "Square Grids" >> $terminal_out
	for file in `ls -v graphs/grid_graphs/*.npz`; do
		echo "Running $file on $node nodes..."
		(time mpiexec -N $node python3 steps.py $file "results/benchmark_grids_$node.csv") &>> $terminal_out
	done

	echo "Random Graphs" >> $terminal_out
	for file in `ls -v graphs/random_graphs/*.npz`; do
		echo "Running $file on $node nodes..."
		(time mpiexec -N $node python3 steps.py $file "results/benchmark_randoms_$node.csv") &>> $terminal_out
	done

	echo "Complete Graphs" >> $terminal_out
	for file in `ls -v graphs/complete_graphs/*.npz`; do
		echo "Running $file on $node nodes..."
		(time mpiexec -N $node python3 steps.py $file "results/benchmark_completes_$node.csv") &>> $terminal_out
	done

done

echo "Benchmarks complete."

python3 plot.py

echo "Plotting complete."
