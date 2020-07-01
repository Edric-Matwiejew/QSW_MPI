#!/bin/bash

python3 cray_graph_gen.py

bash slurm_gen.sh 2 12 04:00:00 "hybrid" "base.slurm" "QSW_MPI_hybrid"
bash slurm_gen.sh 12 1 06:00:00 "pure" "base.slurm" "QSW_MPI_pure"

( cd hybrid ; bash launch.sh )
( cd pure ; bash launch.sh )
