#!/bin/bash

TOP="../../"
BINDIR=$TOP"build/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="hkc_build_pop"
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1x

$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Vaccine_Probability_Susceptible 1  1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a 
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b 



#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 10 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 10 Vaccine_Probability_Susceptible 1  0

#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b Maximum_Number_CT_Vaccinations_Each_Day 1 200 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 30 Vaccine_Probability_Susceptible 1  1

#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1c Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 1 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1d Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  1  Start_Time_Self_Isolation 1 1 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2a Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2b Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2c Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2d Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3a Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3b Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3c Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3d Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 99 Vaccine_Probability_Susceptible  1  1

#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a1 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b1 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1c1 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1d1 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2a1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2b1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2c1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2d1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3a1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3b1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3c1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3d1 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.1 Probability_Trace_Neighbour  1 0.1 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1

#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1a2 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1b2 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1c2 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1d2 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2a2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2b2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2c2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2d2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  0
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3a2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  999  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3b2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  10  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3c2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  15  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3d2 Maximum_Number_CT_Vaccinations_Each_Day 1 30 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour  1 0.5 Start_Time_Contact_Tracing 1  20  Start_Time_Self_Isolation 1 999 Vaccine_Probability_Susceptible 1  1



