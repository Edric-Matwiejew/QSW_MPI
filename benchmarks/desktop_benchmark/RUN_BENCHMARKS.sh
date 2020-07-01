#!/bin/bash

if [ ! -e graphs ]; then
	mkdir graphs
	python3 graphgen.py
fi

export total_time=300

( cd local_step ; bash local_step_benchmark_run.sh )
( cd global_step ; bash global_step_benchmark_run.sh )
( cd series ; bash series_benchmark_run.sh )
