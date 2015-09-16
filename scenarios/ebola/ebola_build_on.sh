#!/bin/bash

TOP="../../"
BINDIR=$TOP"g++/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="ebola_pop"

$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM}

