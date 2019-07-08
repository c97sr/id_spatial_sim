#!/bin/bash

TOP="../../"
BINDIR=$TOP"build/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output/"
NETDIR="./network/"
OUTSTEM="hkc_build_pop"
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1 Maximum_Number_CT_Vaccinations_Each_Day 1 20 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  30  Start_Time_Self_Isolation 1 30 Vaccine_Probability_Susceptible 1  0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2 Maximum_Number_CT_Vaccinations_Each_Day 1 20 Probability_Trace_Household 1 1 Probability_Trace_Neighbour  1 1 Start_Time_Contact_Tracing 1  30  Start_Time_Self_Isolation 1 30 Vaccine_Probability_Susceptible 1  1
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3 Maximum_Number_CT_Vaccinations_Each_Day 1 20 Probability_Trace_Household 1 0 Probability_Trace_Neighbour  1 0 Start_Time_Contact_Tracing 1  30  Start_Time_Self_Isolation 1 30 Vaccine_Probability_Susceptible 1  0
$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_4 Maximum_Number_CT_Vaccinations_Each_Day 1 20 Probability_Trace_Household 1 0 Probability_Trace_Neighbour  1 0 Start_Time_Contact_Tracing 1  30  Start_Time_Self_Isolation 1 30 Vaccine_Probability_Susceptible 1  1



#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_1 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour 1 1 intNoNodes 1 100 dblAverageHousehold  1  3 Vaccine_Probability_Susceptible 1 0 
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_2 Maximum_Number_CT_Vaccinations_Each_Day 1 20000 Probability_Trace_Household 1 1 Probability_Trace_Neighbour 1 1 intNoNodes 1 100 dblAverageHousehold  1  3 Vaccine_Probability_Susceptible  1   1
#$BINDIR$RUN both_run.in $NETDIR${OUTSTEM} ${OUTDIR}ebola_3 Maximum_Number_CT_Vaccinations_Each_Day 1 0 Probability_Trace_Household 1 0.5 Probability_Trace_Neighbour 1 0.5 intNoNodes 1 100 dblAverageHousehold  1  3 Vaccine_Probability_Susceptible  1   1



