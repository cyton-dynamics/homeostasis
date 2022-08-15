#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=50GB

nproc
julia -t `nproc` ./run.jl

