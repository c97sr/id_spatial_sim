#!/bin/bash

TOP="../../"
BINDIR=$TOP"build/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="hkc_build_pop"

$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola


