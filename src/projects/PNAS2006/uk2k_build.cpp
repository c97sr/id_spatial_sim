#include"uk2k.h"
#include<ctime>
// /G7 /O3 /QxW
// /D_WIN32_WINNT=0x0502 /D SR_DEBUG /D SR_BYTEPACKED /D SR_PAE_PAGING

int main(int argc, char* argv[]) {
	
	cerr << sizeof(SR::Node) << "\n\n"; 

	// Read from command line
	int intNoArgs = 2;
	if (argc<intNoArgs+1) SR::srerror("First two arguments parameter_file_name and output_file_name. Rest parsed as parameter values.\n");
	string strParamFile, strOutputFile, strHouseholdFile, strWorkplaceFile, strArgs;
	strParamFile = argv[1];
	strOutputFile = argv[2]; 
	if ((argc-(intNoArgs+1))%3!=0) SR::srerror("An even number of parameter arguments are required.\n");
	for (int i=intNoArgs+1;i<argc;i=i+3){strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";};

	// Set up parameter object
	SR::ParameterSet ukPars;
	cerr  <<  "Reading parameters from file...";
	ukPars.ReadParamsFromFile(strParamFile);
	cerr  <<  "done.\n";
	cerr  <<  "Command line parameters listed below...\n";
	cerr  <<  strArgs << " done.\n";
	if (argc > 4) ukPars.ReadParams(strArgs);
	SetUpExtraParameters(ukPars);
	strHouseholdFile = ukPars.GetTag("strHouseholdDensityFile");
	strWorkplaceFile = ukPars.GetTag("strWorkplaceDensityFile");
	double dblDistanceAllWorkplacesGridDx = ukPars.GetValue("dblDistanceAllWorkplacesGridDx");
	double dblDistanceAllWorkplacesHistDx = ukPars.GetValue("dblDistanceAllWorkplacesHistDx");

	cerr  << "X range from " << ukPars.GetValue("dblXGridMin") << " to " << ukPars.GetValue("dblXGridMin")+ukPars.GetValue("dblXGridSize") << "\n"
		  << "Y range from " << ukPars.GetValue("dblYGridMin") << " to " << ukPars.GetValue("dblYGridMin")+ukPars.GetValue("dblYGridSize") << "\n";

	// Declare stream objects
	ofstream ofs;
	ifstream ifs;
	
	// Initialise the ukPars.RunP.intSeed 
	int tmpseedcheck = ukPars.intSeed;
	ukPars.intSeed = -1*ukPars.intSeed;

	// Initialise basic hexagon and density structures
	SR::Hexagon tmphex(ukPars);
	SR::DensityField PopulationDensityField(strHouseholdFile, 1, 1, ukPars.GetValue("dblXGridMin"), ukPars.GetValue("dblXGridMin")+ukPars.GetValue("dblXGridSize"), ukPars.GetValue("dblYGridMin"), ukPars.GetValue("dblYGridMin")+ukPars.GetValue("dblYGridSize"));
	SR::DensityField WorkplaceDensityField(strWorkplaceFile, 1, 1, ukPars.GetValue("dblXGridMin"), ukPars.GetValue("dblXGridMin")+ukPars.GetValue("dblXGridSize"), ukPars.GetValue("dblYGridMin"), ukPars.GetValue("dblYGridMin")+ukPars.GetValue("dblYGridSize"));
	// SR::DensityField PopulationDensityField(strHouseholdFile);
	// SR::DensityField WorkplaceDensityField(strWorkplaceFile);

	// Generate a GridHex from the vector of nodes
	SR::GridHex ukGridHex(ukPars,tmphex, PopulationDensityField);
	// Generate a workplace using the gridhex

	SR::Workplaces ukWorkplaces(ukGridHex, ukPars, WorkplaceDensityField);

	int intCurrentMCMC=0, intStableRuns;
	int maxMCMC = ukPars.GetIntValue("intMCMCMaxSamplesInMillions")/3;
	double propStable = ukPars.GetValue("dblMCMCProportionResample");
	ostringstream foss;

	// Update the workplaces the required amount
	// Run until an average of 1 million steps has delta ln like not significantly different to 0
	// Then run for an additional propStable*current updates and output the distributions half way through the extra 
	SR::OpenNullFile(strOutputFile+"_initial_commutes.out",ukWorkplaces.PrintDistributionOfCommutes());
	SR::OpenNullFile(strOutputFile+"_initial_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());
	ukWorkplaces.GenerateDistributionOfPossibleCommutes(ukGridHex,ukPars,dblDistanceAllWorkplacesGridDx,
		dblDistanceAllWorkplacesHistDx,strOutputFile+"_distance_all_wps.out");
	cerr << flush;
	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",0,maxMCMC,intCurrentMCMC);
	
	intStableRuns = static_cast<int>(static_cast<double>(intCurrentMCMC)*propStable/2+1);
	if (intCurrentMCMC==0) intStableRuns=0;
	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",intStableRuns,maxMCMC,intCurrentMCMC);	
	foss << "_commutes_" << intCurrentMCMC << ".out"; SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes()); foss.str("");
	SR::OpenNullFile(strOutputFile+"_halfway_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());
	cerr  <<  "Half way through MCMC updates: writing *halfway_commutes.out\n";
	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",intStableRuns,maxMCMC,intCurrentMCMC);
	foss << "_commutes_" << intCurrentMCMC << ".out"; SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes()); foss.str("");
	SR::OpenNullFile(strOutputFile+"_final_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());

	ukPars.ChangeValue("intMCMCMaxSamplesInMillions",intCurrentMCMC);

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
		cerr  <<  "Writing text workplace output...";
			ukWorkplaces.WriteWorkplacesToFile(strOutputFile+"_workplaces.out");
			ukWorkplaces.WriteCommutesToFile(ukGridHex,strOutputFile+"_commutes.out");
		cerr  <<  "done.\n";
	}

	// Log seed for later
	int seedLog = ukPars.intSeed;

	// Set up the network
	ukPars.intSeed = -seedLog;
	SR::GenerateAllNeighbours(ukGridHex,ukPars,procSpatialNeighbourSetup,kernNeighbourSetup,kernIntroSetup,evCountSpatialNeighbour,ukWorkplaces);

	//Allocate Large memory block
	ukPars.ChangeValue("intNumberOfBlocks",ukGridHex.CalculateSizeOfNetwork()/ukPars.GetIntValue("intBlockSize")*102/100+30);
	cerr  <<  "Number of blocks required = " << ukPars.GetIntValue("intNumberOfBlocks") << "\n";
	cerr  <<  "Constructing pages object...";
	SR::PagesForThings<SR::Node> ukEvmemInitial(ukPars.GetIntValue("intBlockSize"),ukPars.GetIntValue("intNumberOfBlocks"));
	cerr  <<  "done.\nReserving memory...";
	ukGridHex.ReserveMemoryForNetwork(&ukEvmemInitial);
	cerr  <<  "done\n.";
	cerr  <<  "Assigning households...";
	ukGridHex.AssignHouseholds();
	cerr  <<  "done\n.";

	// Return to logged seed
	ukPars.intSeed = -seedLog;
	SR::GenerateAllNeighbours(ukGridHex,ukPars,procSpatialNeighbourAdd,kernNeighbourSetup,kernIntroSetup,evAddToNeighbourList,ukWorkplaces); 

	// Define network characteristics so that they can be checked
	ukPars.ChangeValue("dblAveCalcNeighbours",ukGridHex.CalculateAverageSpatialNeighbours());
	ukPars.ChangeValue("dblAveCalcHousehold",ukGridHex.CalculateAverageHousehold());

	// Write parameters, gridhex and network to a file
	cerr  <<  "Writing binary output...";
	ofs.open((strOutputFile+".hex").c_str(),ios::binary);
	ofs << ukPars;
	ofs << ukGridHex;
	ukEvmemInitial.WriteToBinaryFile(ofs,ukGridHex.FirstNode());
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs.close();
	cerr  <<  "done.\n";

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
		cerr  <<  "Writing text output...";
		ukGridHex.WriteArcsToFile(strOutputFile+"_arcs.out");
		ukGridHex.WriteNodeLocationsAndSizesToFile(strOutputFile+"_nodes.out");
		cerr  <<  "done.\n";
	}

	// Dumping parameters to file
	cerr  <<  "Writing parameters to file...";
	SR::OpenNullFile(strOutputFile+"_params.out",ukPars.WriteParams());
	cerr  <<  "done.\n";

	return 0;

}
