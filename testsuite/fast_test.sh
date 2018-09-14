#!/bin/bash

BINDIR="../g++/"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

$BINDIR$BUILD ./params/fast_test_build_params.in  ./output/pop1
# $BINDIR$RUN ./params/fast_test_run_params.in ./output/pop1 ./output/pop1_test1 
