#!/bin/bash

BINDIR="../g++/"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

# rm ./output/*
# rm ./figs/*
# $BINDIR$BUILD ./params/fast_test_build_params.in ./output/small_pop1
$BINDIR$RUN ./params/fast_test_run_params.in  ./output/small_pop1 ./output/ft_sp
# Rscript -e "require(knitr); knit2pdf('fast_test.Rnw')"
