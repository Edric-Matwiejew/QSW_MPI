#!/bin/bash
singularity remote login --tokenfile sylabs-token
singularity build -r qsw_mpi.sif qsw_mpi.def
