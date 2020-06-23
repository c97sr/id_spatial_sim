
#include"ebola.h"
#include<ctime>
#include<limits>

using namespace std;

extern gsl_rng * glob_rng;

int main(int argc, char* argv[]) {

	// Setup global for rngs
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	glob_rng = gsl_rng_alloc (T);

	// Declare stream objects
	ofstream ofs;
	ifstream ifs;
	ostringstream foss;

	// Read enough from command line to load GridHex and Parameters
	int intNoArgs = 2;
	if (argc<intNoArgs+1)
	{
		SR::srerror("First two arguments parameter_file_name and output_file_name. Rest parsed as parameter values.\n");
	}

	string strParamFile, strOutputFile, strWorkplaceFile, strArgs;
	strParamFile = argv[1];
	strOutputFile = argv[2];

	//Load Grid
	ifs.open((strOutputFile+"_initial_gridhex.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary gridhex file.");
	SR::GridHex * ukGridHex;
	ukGridHex = new SR::GridHex();
	ifs >> *ukGridHex;
	ifs.close();

	//Load Parameters
	ifs.open((strOutputFile+"_initial_params.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary parameter file.");
	SR::ParameterSet ukPars;
	ifs >> ukPars;
	ifs.close();

	//Load network generation parameters into ukPars
	if ((argc-(intNoArgs+1))%3!=0)
	{
		SR::srerror("An even number of parameter arguments are required.\n");
	}

	for (int i=intNoArgs+1;i<argc;i=i+3)
	{
		strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";
	}

	ukPars.ReadParamsFromFile(strParamFile);
	if (argc > 4)
	{
		ukPars.ReadParams(strArgs);
	}

	// Mask the nodes in GridHex if needed
	SR::NodeMask mask1(*ukGridHex);
	//int intAgeLB = 10;//ukPars.GetIntValue("intAgeLowerBound");
	//int intAgeUB = 50;//ukPars.GetIntValue("intAgeUpperBound");
	//mask1.AgeMask(intAgeLB,intAgeUB,*ukGridHex);
	mask1.NullMask(*ukGridHex);

	// Set up the workplace density file
	strWorkplaceFile = ukPars.GetTag("strWorkplaceDensityFile");
	double dblDistanceAllWorkplacesGridDx = ukPars.GetValue("dblDistanceAllWorkplacesGridDx");
	double dblDistanceAllWorkplacesHistDx = ukPars.GetValue("dblDistanceAllWorkplacesHistDx");
	SR::DensityField WorkplaceDensityField(strWorkplaceFile);

	SR::Workplaces ukWorkplaces(*ukGridHex, ukPars, WorkplaceDensityField);

	// Assign regional location to nodes and hexagons
	setGZRegionMembership(*ukGridHex);

	int intCurrentMCMC=0, intStableRuns;
	int maxMCMC = ukPars.GetIntValue("intMCMCMaxSamplesInMillions")/3;
	double propStable = ukPars.GetValue("dblMCMCProportionResample");

	// Update the workplaces the required amount
	SR::OpenNullFile(strOutputFile+"_commute_dist_0.csv",ukWorkplaces.PrintDistributionOfCommutes());
	SR::OpenNullFile(strOutputFile+"_initial_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());
	ukWorkplaces.GenerateDistributionOfPossibleCommutes(*ukGridHex,ukPars,dblDistanceAllWorkplacesGridDx,
		dblDistanceAllWorkplacesHistDx,strOutputFile+"_distance_all_wps.out");

	ukWorkplaces.MCMCUpdate(
			*ukGridHex,
			ukPars,
			fnOffsetPower,
			strOutputFile+"_target_commutes_function.out",
			mask1,
			0,
			maxMCMC,
			intCurrentMCMC);

	intStableRuns = static_cast<int>(static_cast<double>(intCurrentMCMC)*propStable/2+1);
	if (intCurrentMCMC==0)
	{
		intStableRuns=0;
	}

	ukWorkplaces.MCMCUpdate(
			*ukGridHex,
			ukPars,
			fnOffsetPower,
			strOutputFile+"_target_commutes_function.out",
			mask1,
			intStableRuns,
			maxMCMC,
			intCurrentMCMC);

	foss << "_commute_dist_" << intCurrentMCMC << ".csv";
	SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes());
	foss.str("");
	SR::OpenNullFile(strOutputFile+"_halfway_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());

	ukWorkplaces.MCMCUpdate(
			*ukGridHex,
			ukPars,
			fnOffsetPower,
			strOutputFile+"_target_commutes_function.out",
			mask1,
			intStableRuns,
			maxMCMC,
			intCurrentMCMC);

	foss << "_commute_dist_" << intCurrentMCMC << ".csv";
	SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes());
	foss.str("");
	SR::OpenNullFile(strOutputFile+"_final_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());

	ukPars.ChangeValue("intMCMCMaxSamplesInMillions",intCurrentMCMC);

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE")
	{
			ukWorkplaces.WriteWorkplacesToFile(strOutputFile+"_workplaces.out");
			ukWorkplaces.WriteCommutesToFile(*ukGridHex,strOutputFile+"_commutes.out");
	}

	// Log seed for later
	int seedLog = ukPars.intSeed;
	unsigned long int checkpoint = gsl_rng_get(glob_rng);
	gsl_rng_set(glob_rng,checkpoint);

	// Set up the network
	ukPars.intSeed = -seedLog;
	SR::GenerateAllNeighbours(*ukGridHex,ukPars,procSpatialNeighbourSetup,kernNeighbourSetup,kernIntroSetup,evCountSpatialNeighbour,ukWorkplaces);


	//Allocate Large memory block
	ukPars.ChangeValue("intNumberOfBlocks",ukGridHex->CalculateSizeOfNetwork()/ukPars.GetIntValue("intBlockSize")*102/100+30);
	SR::PagesForThings<SR::Node> ukEvmemInitial(ukPars.GetIntValue("intBlockSize"),ukPars.GetIntValue("intNumberOfBlocks"));
	ukGridHex->ReserveMemoryForNetwork(&ukEvmemInitial);
	ukGridHex->AssignHouseholds();


	// Return to logged seed
	ukPars.intSeed = -seedLog;
	gsl_rng_set(glob_rng,checkpoint);
	SR::GenerateAllNeighbours(*ukGridHex,ukPars,procSpatialNeighbourAdd,kernNeighbourSetup,kernIntroSetup,evAddToNeighbourList,ukWorkplaces);

	// Define network characteristics so that they can be checked
	ukPars.ChangeValue("dblAveCalcNeighbours",ukGridHex->CalculateAverageSpatialNeighbours());
	ukPars.ChangeValue("dblAveCalcHousehold",ukGridHex->CalculateAverageHousehold());

	// Write parameters, gridhex and network to a file
	ofs.open((strOutputFile+"_params.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs << ukPars;
	ofs.close();
	ofs.open((strOutputFile+"_pages.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ukEvmemInitial.WriteToBinaryFile(ofs,ukGridHex->FirstNode());
	ofs.close();
	ofs.open((strOutputFile+"_gridhex.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs << *ukGridHex;
	ofs.close();


	// Write population density actually used to a file
	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE")
	{
		ukGridHex->WriteArcsToFile(strOutputFile+"_arcs.csv");
		ukGridHex->WriteNodeLocationsAndSizesToFile(strOutputFile+"_nodes.csv");
		cerr  <<  "done.\n";
	}
	// Dumping parameters to file
	SR::OpenNullFile(strOutputFile+"_params.out",ukPars.WriteParams());

	cerr << "Finished.\n";

	// None managed garbage collection
	gsl_rng_free (glob_rng);
	ukGridHex->~GridHex();


	return 0;

}
