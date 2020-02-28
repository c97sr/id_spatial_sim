#!/bin/bash

TOP="../../"
BINDIR=$TOP"g++/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="ebola_pop_monrovia"

$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_1 R0_Network 1	2.71	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_2  R0_Network 1	2.72	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_3  R0_Network 1	2.73	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_4  R0_Network 1	2.74	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_5  R0_Network 1	2.745	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_6  R0_Network 1	2.75	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"} ${OUTDIR}ebola_alpha0_7  R0_Network 1	2.76	0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_8 R0_Network 1	2.67	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_9  R0_Network 1	2.68	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_10  R0_Network 1	2.69	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_11 R0_Network 1	2.7	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_12 R0_Network 1	2.71	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_13 R0_Network 1	2.72	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_14 R0_Network 1	2.73	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_15 R0_Network 1	2.74	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_16 R0_Network 1	2.75	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"} ${OUTDIR}ebola_alpha1_17 R0_Network 1	2.76	1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_18 R0_Network 1	2.71	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_19 R0_Network 1	2.72	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_20 R0_Network 1	2.73	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_21 R0_Network 1	2.74	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_22 R0_Network 1	2.75	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"} ${OUTDIR}ebola_alpha2_23 R0_Network 1	2.76	2
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_24 R0_Network 1	2.71	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_25 R0_Network 1	2.2	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_26 R0_Network 1	2.73	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_27 R0_Network 1	2.74	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_28 R0_Network 1	2.75	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_29 R0_Network 1	2.76	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_30 R0_Network 1	2.77	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_31 R0_Network 1	2.78	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"} ${OUTDIR}ebola_alpha3_32 R0_Network 1	2.79	3
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_33 R0_Network 1	2.82	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_34 R0_Network 1	2.83	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_35 R0_Network 1	2.84	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_36 R0_Network 1	2.85	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_37 R0_Network 1	2.86	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"} ${OUTDIR}ebola_alpha4_38 R0_Network 1	2.87	4
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_39 R0_Network 1	2.84	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_40 R0_Network 1	2.85	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_41 R0_Network 1	2.86	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_42 R0_Network 1	2.87	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_43 R0_Network 1	2.88	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_44 R0_Network 1	2.89	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_45 R0_Network 1	2.9	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_46 R0_Network 1	2.91	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_47 R0_Network 1	2.92	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_48 R0_Network 1	2.93	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"} ${OUTDIR}ebola_alpha5_49 R0_Network 1	2.94	5
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_50 R0_Network 1	2.95	6
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_51 R0_Network 1	2.96	6
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_52 R0_Network 1	2.97	6
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_53 R0_Network 1	2.98	6
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_54 R0_Network 1	2.99	6
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"} ${OUTDIR}ebola_alpha6_55 R0_Network 1	3	6
