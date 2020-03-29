#!/bin/bash

TOP="../../"
BINDIR=$TOP"g++/"
PARSTEM="both"
BUILD="ebola_build.exe"
RUN="ebola_run.exe"

OUTDIR="./output_MCMC/"
NETDIR="./network/"
OUTSTEM="ebola_pop_"

$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha0"}  dblAverageHousehold   1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14										
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha1"}  dblAverageHousehold    1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	1	Commute_Power_One_HK 1	1	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha2"}  dblAverageHousehold    1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	2	Commute_Power_One_HK 1	2	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha3"}  dblAverageHousehold    1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	3	Commute_Power_One_HK 1	3	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha4"}  dblAverageHousehold    1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	4	Commute_Power_One_HK 1	4	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha5"}  dblAverageHousehold   1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	5	Commute_Power_One_HK 1	5	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
$BINDIR$BUILD ${PARSTEM}_setup.in $NETDIR${OUTSTEM+"ebola_pop_alpha6"}  dblAverageHousehold   1	4	dblAverageWorkplaceSize 1	50	dblPropColleguesInNetwork 1	0.14	intMCMCMaxSamplesInMillions 1	210	Commute_Power_One_GZ 1	6	Commute_Power_One_HK 1	6	Commute_Change_Point_GZ 1	0.1	Commute_Change_Point_HK 1	0.1
