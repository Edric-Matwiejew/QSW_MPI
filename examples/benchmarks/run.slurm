#!/bin/bash -l
#SBATCH --job-name=FreeQSW
#SBATCH --partition=workq
#SBATCH --nodes=5
#SBATCH --time=03:00:00
#SBATCH --account=pawsey0309
#SBATCH --output=magnus_bench.out
#SBATCH --export=All

module swap PrgEnv-cray PrgEnv-gnu
module load mpi4py
module load numpy
module load scipy
module load h5py

export OMP_NUM_THREADS=1

mkdir graphs/line_graphs/sym
mkdir graphs/grid_graphs/sym
mkdir graphs/random_graphs/sym
mkdir graphs/complete_graphs/sym

python graphgen.py

echo "FreeQSW benchmarks" > benchmarks.dat

echo "Line Graphs" >> benchmarks.dat
for file in graphs/line_graphs/*.npz
do
	(time srun -N 5 -n 120 -c 1 python3 benchmark_steps.py $file 'benchmark_lines.csv') &>> benchmarks.dat
done

echo "Square Grids" >> benchmarks.dat
for file in graphs/grid_graphs/*.npz
do
	(time srun -N 5 -n 120 -c 1 python3 benchmark_steps.py $file 'benchmark_grids.csv') &>> benchmarks.dat
done

echo "Random Graphs" >> benchmarks.dat
for file in graphs/random_graphs/*.npz
do
	(time srun -N 5 -n 120 -c 1 python3 benchmark_steps.py $file 'benchmark_randoms.csv') &>> benchmarks.dat
done

echo "Complete Graphs" >> benchmarks.dat
for file in graphs/complete_graphs/*.npz
do
	(time srun -N 5 -n 120 -c 1 python3 benchmark_steps.py $file 'benchmark_completes.csv') &>> benchmarks.dat
done
