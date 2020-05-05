#!/bin/bash

TOP="../../"
BINDIR=$TOP"build/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="ebola_pop"

mkdir $OUTDIR

$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM}
