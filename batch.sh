#!/bin/bash

function sighandler() {
  echo
  echo got $1 at `date`
}

trap "sighandler INT" INT
trap "sighandler CONT" CONT
trap "sighandler TERM" TERM

NPROC=`nproc`

echo julia --threads `nproc` ./run.jl --trial $1 --threshold $2 --initial-population $3
julia --threads `nproc` ./run.jl --trial $1 --threshold $2 --initial-population $3

echo Run complete, exit code: $?


