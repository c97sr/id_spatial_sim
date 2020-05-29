#!/bin/bash

BINDIR="../build/"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

mkdir output
echo building...
$BINDIR$BUILD ./params/fast_test_build_params.in  ./output/pop1
$BINDIR$RUN ./params/fast_test_run_params.in ./output/pop1 ./output/pop1_test1 
