#!/bin/bash -l
#SBATCH --job-name=NAME_PLACEHOLDER
#SBATCH --partition=workq
#SBATCH --nodes=NODES_PLACEHOLDER
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --account=
#SBATCH --output=OUTPUT_PLACEHOLDER

export OMP_NUM_THREADS=THREADS_PLACEHOLDER

module load singularity
echo "Starting line graph..."
command time -v srun --export=all -N NODES_PLACEHOLDER -n MPI_PROCS_PLACEHOLDER --cpu_bind=sockets -c THREADS_PLACEHOLDER singularity exec ../qsw_mpi.sif python3 ../cray_QSW_MPI_step.py "../graphs/line_graph_5050.mtx"  "../graphs/sym/line_graph_5050_sym.mtx" RESULTS_PLACEHOLDER_LINE CSV_PLACEHOLDER
echo "...done."
echo "Starting grid graph..."
command time -v srun --export=all -N NODES_PLACEHOLDER -n MPI_PROCS_PLACEHOLDER --cpu_bind=sockets -c THREADS_PLACEHOLDER singularity exec ../qsw_mpi.sif python3 ../cray_QSW_MPI_step.py "../graphs/grid_graph_62.mtx"  "../graphs/sym/grid_graph_62_sym.mtx" RESULTS_PLACEHOLDER_GRID CSV_PLACEHOLDER
echo "...done."
echo "Starting random graph..."
command time -v srun --export=all -N NODES_PLACEHOLDER -n MPI_PROCS_PLACEHOLDER --cpu_bind=sockets -c THREADS_PLACEHOLDER singularity exec ../qsw_mpi.sif python3 ../cray_QSW_MPI_step.py "../graphs/random_graph_2020.mtx"  "../graphs/sym/random_graph_2020_sym.mtx" RESULTS_PLACEHOLDER_RANDOM CSV_PLACEHOLDER
echo "...done."
echo "Starting complete graph..."
command time -v srun --export=all -N NODES_PLACEHOLDER -n MPI_PROCS_PLACEHOLDER --cpu_bind=sockets -c THREADS_PLACEHOLDER singularity exec ../qsw_mpi.sif python3 ../cray_QSW_MPI_step.py "../graphs/complete_graph_400.mtx"  "../graphs/sym/complete_graph_400_sym.mtx" RESULTS_PLACEHOLDER_COMPLETE CSV_PLACEHOLDER
echo "...done."
