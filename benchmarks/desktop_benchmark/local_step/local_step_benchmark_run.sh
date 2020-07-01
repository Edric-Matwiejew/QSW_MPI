#!/bin/bash

if [ -e Results ]; then
	rm -r Results
	rm -r Plots
fi

mkdir -p Results
mkdir -p Plots

export OMP_NUM_THREADS=1
nodes=(1 4)

line_args="$* "../graphs/line_graphs/" "../graphs/line_graphs/sym/" "Results/line_graphs" $total_time"
grid_args="$* "../graphs/grid_graphs/" "../graphs/grid_graphs/sym/" "Results/grid_graphs" $total_time"
random_args="$* "../graphs/random_graphs/" "../graphs/random_graphs/sym/" "Results/random_graphs" $total_time"
complete_args="$* "../graphs/complete_graphs/" "../graphs/complete_graphs/sym/" "Results/complete_graphs" $total_time"

echo "Started QSWalk.jl local steps..."
julia QSWalk_jl_local_step.jl $line_args
julia QSWalk_jl_local_step.jl $grid_args
julia QSWalk_jl_local_step.jl $random_args
julia QSWalk_jl_local_step.jl $complete_args
echo "Finished."

echo "Started QSWalk.m local steps..."
wolframscript -f QSWalk_m_local_step.wls $line_args
wolframscript -f QSWalk_m_local_step.wls $grid_args
wolframscript -f QSWalk_m_local_step.wls $random_args
wolframscript -f QSWalk_m_local_step.wls $complete_args
echo "Finished."

echo "Started QSW_MPI local steps..."
for node in ${nodes[@]}; do
	mpiexec -N $node python3 QSW_MPI_local_step.py $line_args
	mpiexec -N $node python3 QSW_MPI_local_step.py $grid_args
	mpiexec -N $node python3 QSW_MPI_local_step.py $random_args
	mpiexec -N $node python3 QSW_MPI_local_step.py $complete_args
	echo "$node done."
done

python3 plot_accuracy.py
python3 plot_times.py
