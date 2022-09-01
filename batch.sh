#!/bin/bash

NPROC=`nproc`

echo julia -t `nproc` ./run.jl --trial $1 --threshold $2 --initial-population $3
julia -t `nproc` ./run.jl --trial $1 --threshold $2 --initial-population $3

echo Run complete, exit code: $?


