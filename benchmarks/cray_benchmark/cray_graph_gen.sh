#!/bin/bash

# To be ran from an interactive session.

module load singularity

singularity exec qsw_mpi.sif python3 cray_graph_gen.py


