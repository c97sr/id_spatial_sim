#!/bin/bash

BINDIR="../g++/"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

$BINDIR$BUILD ./pop1_params.in ./output/pop1
$BINDIR$RUN ./output/pop1 ./output/pop1_test1 
