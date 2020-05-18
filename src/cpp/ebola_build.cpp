/*  Copyright 2018 Steven Riley.

    This file is part of id_spatial_sim.

    id_spatial_sim is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    id_spatial_sim is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with id_spatial_sim.  If not, see <https://www.gnu.org/licenses/>.  */

// Just testing the web interface for git hub
// And also the add and commit functions

#include"ebola.h"
#include<ctime>
#include<limits>

using namespace std;

// Effectively a global for the random number generator
// There is a job to find all the uses of this and then add explicit arguments for the random number
// object to be passed up and down the call stack. But for now, its not the end of the world to
// have it as a global
extern gsl_rng * glob_rng;

int main(int argc, char* argv[]) {

	// Setup the global random number generator
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	glob_rng = gsl_rng_alloc (T);

	// Declare stream objects
	ofstream ofs;
	ifstream ifs;
	ostringstream foss;

	// Read from command line and set up parameter object
	int intNoArgs = 2;
	if (argc<intNoArgs+1) SR::srerror("First two arguments parameter_file_name and output_file_name. Rest parsed as parameter values.\n");
	string strParamFile, strOutputFile, strHouseholdFile, strWorkplaceFile, strArgs;
	strParamFile = argv[1];
	strOutputFile = argv[2];
	if ((argc-(intNoArgs+1))%3!=0) SR::srerror("An even number of parameter arguments are required.\n");
	for (int i=intNoArgs+1;i<argc;i=i+3){strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";};
	SR::ParameterSet ukPars;
	ukPars.ReadParamsFromFile(strParamFile);
	if (argc > 4) ukPars.ReadParams(strArgs);

	// Set up the population and workplace density files
	strHouseholdFile = ukPars.GetTag("strHouseholdDensityFile");
	strWorkplaceFile = ukPars.GetTag("strWorkplaceDensityFile");
	SR::DensityField PopulationDensityField(strHouseholdFile);
	// PopulationDensityField.WriteAsciiGrid("tmpAsciiDump.out");
	SR::DensityField WorkplaceDensityField(strWorkplaceFile);
	ukPars.AddValue("dblXGridMin",PopulationDensityField.GetMinX());
	ukPars.AddValue("dblXGridSize",PopulationDensityField.GetMaxX()-PopulationDensityField.GetMinX());
	ukPars.AddValue("dblYGridMin",PopulationDensityField.GetMinY());
	ukPars.AddValue("dblYGridSize",PopulationDensityField.GetMaxY()-PopulationDensityField.GetMinY());

	SetUpExtraParameters(ukPars);

	double dblDistanceAllWorkplacesGridDx = ukPars.GetValue("dblDistanceAllWorkplacesGridDx");
	double dblDistanceAllWorkplacesHistDx = ukPars.GetValue("dblDistanceAllWorkplacesHistDx");

	// Initialise the ukPars.RunP.intSeed
	ukPars.intSeed = -1*ukPars.intSeed;
	gsl_rng_set(glob_rng, ukPars.intSeed);

	// Initialise basic hexagon and density structures
	SR::Hexagon tmphex(ukPars);

	// Generate a GridHex from the vector of nodes
	// Error here. Need a pointer to make this work, which may be a pain!
	SR::GridHex * ukGridHex;

	string strExistingNetworkFile = ukPars.GetTag("fnExistingNetwork");
	if (strExistingNetworkFile == "FALSE") {

		ukGridHex = new SR::GridHex(ukPars,tmphex, PopulationDensityField);

	} else {

		// Load gridhex from a binary
		cerr << "Opening binary file for existing gridhex and checking for consistency...\n";
		ifs.open((strExistingNetworkFile).c_str(),ios::binary);
		if (ifs.fail()) SR::srerror("Problem opening binary gridhex file.");
		ukGridHex = new SR::GridHex(ukPars,tmphex,ifs);;
		ifs.close();

	}

	// Mask the nodes in GridHex if needed
	SR::NodeMask mask1(*ukGridHex);
	int intAgeLB = ukPars.GetIntValue("intAgeLowerBound");
	int intAgeUB = ukPars.GetIntValue("intAgeUpperBound");
	mask1.AgeMask(intAgeLB,intAgeUB,*ukGridHex);
	// mask1.NullMask(*ukGridHex);

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
	if (intCurrentMCMC==0) intStableRuns=0;

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

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
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

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
		ukGridHex->WriteArcsToFile(strOutputFile+"_arcs.csv");
		ukGridHex->WriteNodeLocationsAndSizesToFile(strOutputFile+"_nodes.csv");
		cerr  <<  "done.\n";
	}

	// Dumping parameters to file
	SR::OpenNullFile(strOutputFile+"_params.out",ukPars.WriteParams());

	cerr << "Finished.\n";

	// Now debug routine to make GridHex roundtrip safe


	// None managed garbage collection
	gsl_rng_free (glob_rng);
	ukGridHex->~GridHex();

	return 0;

}
