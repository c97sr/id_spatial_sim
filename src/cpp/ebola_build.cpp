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

extern gsl_rng * glob_rng;


// Some globals ONLY ONLY for random numbers
// I hate globals, but they seem to make sense here
// gsl_rng * glob_rng;

// Probably need an object for multiple travel behaviors that can handle many different zones
// First thing is a multiple cached dbl object

// This function needs to recreate a numerical recepies

int main(int argc, char* argv[]) {

    // Example code for the setup of a global GSL random number generator
	// Setup random number generator before ever used
	// The pattern of externals in the linked files is not obvious to me
	// It needs to be an external here and then a simple global in all the otehr files
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	glob_rng = gsl_rng_alloc (T);

	gsl_rng_set(glob_rng, 123456);
	cout << SR::rngtest() << endl;
	cout << endl << endl;
	unsigned long int check = gsl_rng_get(glob_rng);
	gsl_rng_set(glob_rng,check);
    cout << SR::rngtest() << endl;
    cout << endl << endl;
	gsl_rng_set(glob_rng, check);
    cout << SR::rngtest() << endl;

    // return 0;

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

	// Declare stream objects
        ofstream ofs;
	ifstream ifs;
	ostringstream foss;

	// Initialise the ukPars.RunP.intSeed
	ukPars.intSeed = -1*ukPars.intSeed;
	gsl_rng_set(glob_rng, ukPars.intSeed);

	// Initialise basic hexagon and density structures
	SR::Hexagon tmphex(ukPars);

	// Generate a GridHex from the vector of nodes
	SR::GridHex ukGridHex(ukPars,tmphex, PopulationDensityField);
	SR::Workplaces ukWorkplaces(ukGridHex, ukPars, WorkplaceDensityField);

	// Assign regional location to nodes and hexagons
	setGZRegionMembership(ukGridHex);

	int intCurrentMCMC=0, intStableRuns;
	int maxMCMC = ukPars.GetIntValue("intMCMCMaxSamplesInMillions")/3;
	double propStable = ukPars.GetValue("dblMCMCProportionResample");

	// Update the workplaces the required amount
	SR::OpenNullFile(strOutputFile+"_initial_commutes.out",ukWorkplaces.PrintDistributionOfCommutes());
	SR::OpenNullFile(strOutputFile+"_initial_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());
	ukWorkplaces.GenerateDistributionOfPossibleCommutes(ukGridHex,ukPars,dblDistanceAllWorkplacesGridDx,
		dblDistanceAllWorkplacesHistDx,strOutputFile+"_distance_all_wps.out");

	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",0,maxMCMC,intCurrentMCMC);

	intStableRuns = static_cast<int>(static_cast<double>(intCurrentMCMC)*propStable/2+1);
	if (intCurrentMCMC==0) intStableRuns=0;

	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",intStableRuns,maxMCMC,intCurrentMCMC);
	foss << "_commutes_" << intCurrentMCMC << ".out";
	SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes());
	foss.str("");
	SR::OpenNullFile(strOutputFile+"_halfway_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());
	ukWorkplaces.MCMCUpdate(ukGridHex, ukPars, fnOffsetPower, strOutputFile+"_target_commutes_function.out",intStableRuns,maxMCMC,intCurrentMCMC);
	foss << "_commutes_" << intCurrentMCMC << ".out";
	SR::OpenNullFile(strOutputFile+foss.str(),ukWorkplaces.PrintDistributionOfCommutes());
	foss.str("");
	SR::OpenNullFile(strOutputFile+"_final_wp_sizes.out",ukWorkplaces.PrintDistributionOfSizes());

	ukPars.ChangeValue("intMCMCMaxSamplesInMillions",intCurrentMCMC);

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
			ukWorkplaces.WriteWorkplacesToFile(strOutputFile+"_workplaces.out");
			ukWorkplaces.WriteCommutesToFile(ukGridHex,strOutputFile+"_commutes.out");
	}

	// Log seed for later
	int seedLog = ukPars.intSeed;
	unsigned long int checkpoint = gsl_rng_get(glob_rng);
	gsl_rng_set(glob_rng,checkpoint);

	// Set up the network
	ukPars.intSeed = -seedLog;
	SR::GenerateAllNeighbours(ukGridHex,ukPars,procSpatialNeighbourSetup,kernNeighbourSetup,kernIntroSetup,evCountSpatialNeighbour,ukWorkplaces);

	//Allocate Large memory block
	ukPars.ChangeValue("intNumberOfBlocks",ukGridHex.CalculateSizeOfNetwork()/ukPars.GetIntValue("intBlockSize")*102/100+30);
	SR::PagesForThings<SR::Node> ukEvmemInitial(ukPars.GetIntValue("intBlockSize"),ukPars.GetIntValue("intNumberOfBlocks"));
	ukGridHex.ReserveMemoryForNetwork(&ukEvmemInitial);
	ukGridHex.AssignHouseholds();

	// Return to logged seed
	ukPars.intSeed = -seedLog;
	gsl_rng_set(glob_rng,checkpoint);
	SR::GenerateAllNeighbours(ukGridHex,ukPars,procSpatialNeighbourAdd,kernNeighbourSetup,kernIntroSetup,evAddToNeighbourList,ukWorkplaces);

	// Define network characteristics so that they can be checked
	ukPars.ChangeValue("dblAveCalcNeighbours",ukGridHex.CalculateAverageSpatialNeighbours());
	ukPars.ChangeValue("dblAveCalcHousehold",ukGridHex.CalculateAverageHousehold());

	// Write parameters, gridhex and network to a file
	ofs.open((strOutputFile+"_params.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror((strOutputFile+"_params.hex - failed to open").c_str());
	ofs << ukPars;
	ofs.close();
	ofs.open((strOutputFile+"_pages.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror((strOutputFile+"_pages.hex - failed to open").c_str());
	ukEvmemInitial.WriteToBinaryFile(ofs,ukGridHex.FirstNode());
	ofs.close();
	ofs.open((strOutputFile+"_gridhex.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror((strOutputFile+"_gridhex.hex - failed to open").c_str());
	ofs << ukGridHex;
	ofs.close();

	// Write population density actually used to a file

	if (ukPars.GetTag("blNetworkDumpFile")=="TRUE") {
		ukGridHex.WriteArcsToFile(strOutputFile+"_arcs.out");
		ukGridHex.WriteNodeLocationsAndSizesToFile(strOutputFile+"_nodes.out");
		cerr  <<  "done.\n";
	}

	// Dumping parameters to file
	SR::OpenNullFile(strOutputFile+"_params.out",ukPars.WriteParams());

	cerr << "Finished.\n";

	gsl_rng_free (glob_rng);

	return 0;

}
