#!/bin/bash

TOP="../../"
BINDIR=$TOP"build/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="hkc_build_pop"

$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM}
