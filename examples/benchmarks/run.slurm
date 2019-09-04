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

echo "1"
time srun --export=ALL -N 5 -n 120 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "2"
time srun --export=ALL -N 5 -n 64 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "4"
time srun --export=ALL -N 5 -n 32 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "8"
time srun --export=ALL -N 5 -n 16 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "16"
time srun --export=ALL -N 5 -n 8 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "32"
time srun --export=ALL -N 5 -n 4 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "64"
time srun --export=ALL -N 5 -n 2 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py
echo "120"
time srun --export=ALL -N 5 -n 1 -c 1 --cpu_bind=sockets mprof run --include-children magnus_bench.py


