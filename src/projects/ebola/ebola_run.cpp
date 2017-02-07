// Compile options to consider
// /G7 /O3 /QxW
// /D_WIN32_WINNT=0x0502 /D SR_BYTEPACKED /D SR_PAE_PAGING

#include"ebola.h"

int main(int argc, char* argv[]) {
  cerr << "\nBuilt on 15/02/2006 at 11:00\n";
#ifdef SR_PAE_PAGING
  cerr << "MS Windows Advanced Server or better required";
  cerr << "to run this version.\n";
#endif

#ifdef SR_BYTEPACKED
  cerr << "Built with a 32 bit byte packed node.\n";
#endif
  
  // Read from command line
  int intNoArgs = 3;
  if (argc<intNoArgs+1) SR::srerror("First three arguments parameter_file_name, binary_file_name, output_file_name. Rest parsed as parameter values.\n");
  string strParamFile, strBinaryFile, strOutputFile, strRunOutputFile, strJunk, strArgs, strValuesChangesFile, strEventFile;
  ostringstream osstmp;
  bool blReset;
  strParamFile = argv[1];
  strBinaryFile = argv[2];
  strOutputFile = argv[3];

	if ((argc-(intNoArgs+1))%3!=0) SR::srerror("An even number of parameter arguments are required.\n");
	for (int i=intNoArgs+1;i<argc;i=i+3){strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";};

	// Define some utility variables
	ofstream ofs;

	// Debug binary in and out
	cerr << "Opening binary file and checking for consistency...\n";
	ifstream ifs;
	ifs.open((strBinaryFile+"_params.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary parameter file.");
	SR::ParameterSet ukPars;
	ifs >> ukPars;
	ifs.close();
	SR::Hexagon tmphex(ukPars.GetIntValue("intNoCharacteristics"),ukPars.GetIntValue("intNoMaximalNodes"));
	ifs.open((strBinaryFile+"_gridhex.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary gridhex file.");
	SR::GridHex ukGridHex(ukPars,tmphex,ifs);
	ifs.close();
	ifs.open((strBinaryFile+"_pages.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary pages file.");
	SR::PagesForThings<SR::Node> ukEvmemInitial(ukPars.GetIntValue("intBlockSizeInUnitsOfPointers"),ukPars.GetIntValue("intNumberOfBlocks"),ifs,ukGridHex.FirstNode());
	ifs.close();
	ukGridHex.AssignToPagesOfNodes(&ukEvmemInitial);

	// cerr << ukGridHex.FirstHexagon()->OutputIndexByCharacteristicOnOneLine() << "\n";

	// Define minimum set of parameters initialize in the absence of parameter file reading
	cerr << "Command line parameters listed below...\n";
	ukPars.ReadParamsFromFile(strParamFile);
	if (argc > 3) ukPars.ReadParams(strArgs);
	RunExtraParameters(ukPars);
	cerr << "done.\n";
	bool blPartialRun;
	if (ukPars.GetTag("strPartialRun")=="TRUE") blPartialRun=true; else blPartialRun=false;
	bool blLogOutputBinary;
	if (ukPars.GetTag("strLogOutputBinary")=="TRUE") blLogOutputBinary=true; else blLogOutputBinary=false;

	// Output text version of network if required
	if (ukPars.GetTag("blNetworkDumpFileRunCode")=="TRUE") {
		cerr << "Writing text output...";
		ukGridHex.WriteArcsToFile(strOutputFile+"_runcode_arcs.out");
		ukGridHex.WriteNodeLocationsAndSizesToFile(strOutputFile+"_runcode_nodes.out");
		cerr << "done.\n";
	}

	// Set up the infection process
	// Characteristic go 0-sus, 1-lat, 2-inf nonsym (prodromal?), 3-inf sym, 4-rec, 5-dead
	// Arg 2 of vector constructor gives value (arg one gives size)
	vector<int> SourceCharacteristic(1,2);
	// vector<int> TargetCharacteristic(1,0);
	vector<int> TargetCharacteristic(9);
	TargetCharacteristic[0]=0;
	TargetCharacteristic[1]=1;
	TargetCharacteristic[2]=2;
	TargetCharacteristic[3]=3;
	TargetCharacteristic[4]=4;
	TargetCharacteristic[5]=5;
	TargetCharacteristic[6]=6;
	TargetCharacteristic[7]=7;
	TargetCharacteristic[8]=8;
	vector<int> SourceCharacteristicAllSymp(2);
	SourceCharacteristicAllSymp[0]=3;
	SourceCharacteristicAllSymp[1]=8;
	vector<int> SourceCharacteristicAllSpatial(2);
	SourceCharacteristicAllSpatial[0]=2;
	SourceCharacteristicAllSpatial[1]=3;
	vector<int> SourceCharacteristicEarlySymp(1,3);
	vector<int> SourceCharacteristicRegionalVaccination(3);
	SourceCharacteristicRegionalVaccination[0]=0;
	SourceCharacteristicRegionalVaccination[1]=1;
	SourceCharacteristicRegionalVaccination[2]=2;
	double *ptConstSpatial = ukPars.GetPointer("Relative_Transmit_Spatial");

	strValuesChangesFile = ukPars.GetTag("strFileParamChanges");
	SR::ParameterValueSetB objParamChanges(strValuesChangesFile);

	if (blPartialRun) strEventFile = ukPars.GetTag("strEventFile");

	// Start parameter loop
	for (int i=0;i<objParamChanges.NoChanges();++i) {

		// ukPars.intSeed = -1234;
		// The line above does not produce the desired results
		// for now assume OK
		// in the future at some point, have an array to log values at different points in the loop and throw an error as soon as it fails

		osstmp.str("");
		osstmp << strOutputFile << "_pset_" << i;
		strRunOutputFile = osstmp.str();

		objParamChanges.UpdateParameterSet(ukPars);

		// Set up monitoring object
		GatherUKPoxInfoA gatherA(static_cast<int>(ukPars.GetValue("dblEndTime")/ukPars.GetValue("dblTimeStep")+1),
			ukPars.GetValue("dblTimeStep"),
			ukPars.GetIntValue("intRealisationsPerParameterSet"),
			strRunOutputFile,ukPars.GetIntValue("intEventLogStackSize"),
			0);

		if (ukPars.GetTag("strEventLogging")=="TRUE") gatherA.StartEventLogging();
		else gatherA.StopEventLogging();

		// Have removed cached kernel from below
		SR::SpatialKernel ukPoxSpatial(
			SourceCharacteristicAllSpatial,
			TargetCharacteristic,
			kernSpatialInfectionOff,
			procInfection,
			ukPars.GetIntValue("intInfectionKernelStackSize"),
			0);
		SR::HouseholdKernel ukPoxHousehold(
			SourceCharacteristic,
			TargetCharacteristic,
			kernHouseholdInfection,
			procInfection,
			0);
		SR::NeighbourKernel ukPoxNeighbour(
			SourceCharacteristic,
			TargetCharacteristic,
			kernNeighbourInfection,
			procInfection,
			0);
		SR::HouseholdKernel ukPoxHouseholdSymp(
			SourceCharacteristicAllSymp,
			TargetCharacteristic,
			kernHouseholdInfectionSymp,
			procInfection,
			0);
		SR::NeighbourKernel ukPoxNeighbourSymp(
			SourceCharacteristicEarlySymp,
			TargetCharacteristic,
			kernNeighbourInfectionSymp,
			procInfection,
			0);
		SR::HexKernel ukRegionalVaccination(
			SourceCharacteristicRegionalVaccination,
			procRegionalVaccination,
			ukGridHex.GetNoHexagons(),
			ukPars.GetIntValue("Max_Daily_Global"),
			ukPars.GetIntValue("Max_Global_Overall"),
			ukPars.GetValue("Max_Local_Per_Person"));

		blReset = static_cast<bool>(ukPars.GetValue("boolResetICEachRealisation"));
		double dblMaxNoInfections = ukPars.GetValue("intMaxCumInfections");

		CalcBetasWithAttackRate(ukPars,ukGridHex,kernSpatialEbola);

		// Declare a suitable event matrix and initialize fixed initial node
		int remMaxDelay=static_cast<int>(ukPars.GetValue("Latent_Vaccinated_One_Maximum_Time")/ukPars.GetValue("dblTimeStep")+1);
		double maxVaccinesPerday=ukPars.GetIntValue("Max_Daily_Global");
		SR::EventMatrix ukEm(ukPars.GetIntValue("intMaxDelayTimesteps"),ukPars.GetIntValue("intMaxDailyEvents"),evVaccinate,static_cast<int>(ukPars.GetValue("Maximum_Number_CT_Vaccinations_Each_Day")));
		SR::EventMatrix ukRegionalEm(remMaxDelay,static_cast<int>(maxVaccinesPerday*1.1),evVaccinate,static_cast<int>(maxVaccinesPerday));

		int intCurrentRealisation,intCurrentTimeStep;
		double dblCurrentCumInfections;
		int intFreqIncRecord = ukPars.GetIntValue("intFrequencyWriteIncidence");
		int intTotalRealisations = ukPars.GetIntValue("intRealisationsPerParameterSet");
		if (intFreqIncRecord > intTotalRealisations) intFreqIncRecord = intTotalRealisations;

		// Set times for regional intervention
		double dblStartRV = ukPars.GetValue("Start_Time_Regional_Vaccination");
		double dblRangeRV = ukPars.GetValue("Range_Regional_Vaccination");
		double dblMaxSympRV = ukPars.GetValue("Max_Symp_Prev");
		double dblMaxSusRV = ukPars.GetValue("Max_Sus_Prev");

		if (ukPars.GetValue("R0_Overall")<1.0 || ukPars.GetValue("Theta")>ukPars.GetValue("dblThetaMax")) intCurrentRealisation = intTotalRealisations;

		int intCountNumberReachingMaxInfections;
		double dblSumTimesMaxInfections;
		int intCountNumberReachingZeroInfectious;
		double dblSumTimesZeroInfectious;
		int intSumInfectionsReachingMaximumTime;
		int intAverageInfectionsIfControlled=0;
		int intAverageInfectionsIfNotControlledMaxNotReached=0;

		int intCurrentlyInfectious;

		intCountNumberReachingMaxInfections=0;
		dblSumTimesMaxInfections=0;
		intCountNumberReachingZeroInfectious=0;
		dblSumTimesZeroInfectious=0;
		intSumInfectionsReachingMaximumTime=0;

		SR::InitialConditions firstIc(ukGridHex,
			ukPars,
			ukPars.GetIntValue("Seeding_number"),
			ukPars.GetIntValue("Max_seeding_number"),
			ukPars.GetValue("Seeding_x"),
			ukPars.GetValue("Seeding_y"),
			ukPars.GetValue("Seeding_range"),
			static_cast<int>(ukPars.GetValue("Seeding_max_try")));

		double dblTrickleRate = ukPars.GetValue("Trickle_rate")/ukPars.GetValue("dblTimeStep");
		double dblTrickleDuration = ukPars.GetValue("Trickle_duration");

		// Start main loop
		for (intCurrentRealisation=0;intCurrentRealisation < intTotalRealisations;++intCurrentRealisation) {
			dblCurrentCumInfections=0;
			intCurrentTimeStep=0;
			ukEm.ClearAllEvents();
			ukPars.ChangeValue("dblCurrentTime",ukPars.GetValue("dblStartTime"));
			ukRegionalVaccination.ResetRegion(ukGridHex);
			if (!blPartialRun) ukGridHex.MakeAllNodesThisCharacteristicAndGeneration(ukPars.GetIntValue("intIndexOfSusceptible"),0);
			if (!blPartialRun) ukGridHex.MakeAllHexagonsNotInRegionalTreatment();
			if (blReset && !blPartialRun) firstIc.Reselect(ukGridHex,ukPars);
			if (!blPartialRun) firstIc.ApplySeed(ukGridHex,ukPars,procInfection,ukEm);
			gatherA.ClearCumPrev();
			gatherA.ResetTimestep();
			if (blPartialRun) ukEm.LoadPendingEventsFromFile(strEventFile,ukGridHex,PoxIntToEvent);
			ukEm.ApplyEvents(ukPars,gatherA);
			gatherA.IncrementTimestep();
			ukPars.ChangeValue("dblCurrentTime", ukPars.GetValue("dblCurrentTime")+ukPars.GetValue("dblTimeStep"));intCurrentTimeStep++;
			while (ukPars.GetValue("dblCurrentTime") <= ukPars.GetValue("dblEndTime")) {
				ukPoxSpatial.GenerateAllEventsWithCross(ukPars,ukEm,ukGridHex,distTestDist,ptConstSpatial);
				ukPoxHousehold.GenerateAllEvents(ukPars,ukEm,ukGridHex);
				ukPoxNeighbour.GenerateAllEvents(ukPars,ukEm,ukGridHex);
				ukPoxHouseholdSymp.GenerateAllEvents(ukPars,ukEm,ukGridHex);
				ukPoxNeighbourSymp.GenerateAllEvents(ukPars,ukEm,ukGridHex);
				ukRegionalVaccination.GenerateEvents(ukRegionalEm,ukPars);
				if (dblTrickleRate > 0) firstIc.Trickle(ukGridHex,ukPars,procInfection,ukEm,dblTrickleRate,dblTrickleDuration);
				ukEm.ApplyEvents(ukPars,gatherA);
				ukEm.StringOutputTheStack();
				ukRegionalEm.ApplyEvents(ukPars,gatherA);
				if (ukPars.GetValue("dblCurrentTime") > dblStartRV)
					ukRegionalVaccination.UpdateRegion(ukGridHex,dblRangeRV,dblMaxSympRV,dblMaxSusRV);
				dblCurrentCumInfections+=*(gatherA.AccessData(0,intCurrentTimeStep,intCurrentRealisation));
				cerr << "R " << intCurrentRealisation << "\tt " << ukPars.GetValue("dblCurrentTime") << "\t"
					 << "In " << *(gatherA.AccessData(0,intCurrentTimeStep,intCurrentRealisation)) << "\t"
					 << "Fe " << *(gatherA.AccessData(4,intCurrentTimeStep,intCurrentRealisation)) << "\t"
					 << "Ra " << *(gatherA.AccessData(5,intCurrentTimeStep,intCurrentRealisation)) << "\t"
					 << "Re " << *(gatherA.AccessData(6,intCurrentTimeStep,intCurrentRealisation)) << "\t"
					 << "Va " << *(gatherA.AccessData(1,intCurrentTimeStep,intCurrentRealisation)) << "     \r";
				if (dblCurrentCumInfections>dblMaxNoInfections) {
					intCountNumberReachingMaxInfections++;
					dblSumTimesMaxInfections+=ukPars.GetValue("dblCurrentTime");
					ukPars.ChangeValue("dblCurrentTime",1e10);
				};
				intCurrentlyInfectious = ukGridHex.GetTotalInfectedNotInfectious();
				if (ukGridHex.GetTotalInfectedNotInfectious()==0 &&
					(dblTrickleRate < 0.000000001 || ukPars.GetValue("dblCurrentTime") > dblTrickleDuration)) {
					intCountNumberReachingZeroInfectious++;
					intAverageInfectionsIfControlled+=static_cast<int>(dblCurrentCumInfections);
					dblSumTimesZeroInfectious+=ukPars.GetValue("dblCurrentTime");
					ukPars.ChangeValue("dblCurrentTime",1e10);
				};
				gatherA.IncrementTimestep(); ukPars.ChangeValue("dblCurrentTime",ukPars.GetValue("dblCurrentTime")+ukPars.GetValue("dblTimeStep"));intCurrentTimeStep++;
			}
			cerr << "\n";
			if (ukPars.GetValue("dblCurrentTime") < 1e10) {
				intSumInfectionsReachingMaximumTime += intCurrentlyInfectious;
				intAverageInfectionsIfNotControlledMaxNotReached += static_cast<int>(dblCurrentCumInfections);
			}
			gatherA.IncrementRealisations();
			if ((intCurrentRealisation+1)%intFreqIncRecord==0)
				if (!gatherA.WriteAllEventIncidencesToFile())
					cerr << "Failed to write incidence file in realisation " << intCurrentRealisation << "\n";
			 }
		if(!gatherA.FlushStack()) {
			cerr << gatherA.GetFileBase() << "_END_\n";
				SR::srerror("Unable to flush stack in GatherUKPoxInfoA::RegisterEventAfterApplication");
		}

		// Dumping parameters to file
		cerr << "Writing parameters to file...";
		ukPars.ChangeValue("dblPropControlled",static_cast<double>(intCountNumberReachingZeroInfectious)/static_cast<double>(intTotalRealisations));
		ukPars.ChangeValue("dblPropReachedMax",static_cast<double>(intCountNumberReachingMaxInfections)/static_cast<double>(intTotalRealisations));
		ukPars.ChangeValue("dblPropNotControlledNotReachedMax",static_cast<double>(intTotalRealisations-intCountNumberReachingZeroInfectious-intCountNumberReachingMaxInfections)/static_cast<double>(intTotalRealisations));
 		if (intCountNumberReachingZeroInfectious>0) {
			ukPars.ChangeValue("dblAveTimeThoseControlled",dblSumTimesZeroInfectious/static_cast<double>(intCountNumberReachingZeroInfectious));
			ukPars.ChangeValue("intAverageInfectionsIfControlled",intAverageInfectionsIfControlled/static_cast<double>(intCountNumberReachingZeroInfectious));
		}
		else {
			ukPars.ChangeValue("dblAveTimeThoseControlled",0);
			ukPars.ChangeValue("intAverageInfectionsIfControlled",0);
		}
		if (intCountNumberReachingMaxInfections>0) {
			ukPars.ChangeValue("dblAveTimeThoseReachedMax",dblSumTimesMaxInfections/static_cast<double>(intCountNumberReachingMaxInfections));
		}
		else {
			ukPars.ChangeValue("dblAveTimeThoseReachedMax",0);
		}
		if (intCountNumberReachingZeroInfectious+intCountNumberReachingMaxInfections != intTotalRealisations) {
			ukPars.ChangeValue("dblAveInfAtEndTimeThoseNotControlledNotReachedMax",
				static_cast<double>(intSumInfectionsReachingMaximumTime)/static_cast<double>(intTotalRealisations-intCountNumberReachingZeroInfectious-intCountNumberReachingMaxInfections));
				ukPars.ChangeValue("intAverageInfectionsIfNotControlledMaxNotReached",
					static_cast<double>(intAverageInfectionsIfNotControlledMaxNotReached)/static_cast<double>(intTotalRealisations-intCountNumberReachingZeroInfectious-intCountNumberReachingMaxInfections));
		} else {
			ukPars.ChangeValue("dblAveInfAtEndTimeThoseNotControlledNotReachedMax",0);
			ukPars.ChangeValue("intAverageInfectionsIfNotControlledMaxNotReached",0);
		}
		SR::OpenNullFile(strRunOutputFile+"_params.out",ukPars.WriteParams());
		ukEm.SaveStackToFile((strRunOutputFile+"_pend_ev.out").c_str(),ukGridHex,PoxEventToInt);
		cerr << "done.\n";
	}

	if (blLogOutputBinary) {
		// Write parameters, gridhex and network to a file
		cerr  <<  "Writing binary output...";
		ofs.open((strOutputFile+"_out.hex").c_str(),ios::binary);
		ofs << ukPars;
		ofs << ukGridHex;
		ukEvmemInitial.WriteToBinaryFile(ofs,ukGridHex.FirstNode());
		if (ofs.fail()) SR::srerror("You idiot.");
		ofs.close();
		cerr  <<  "done.\n";
	};

	return 0;

}
