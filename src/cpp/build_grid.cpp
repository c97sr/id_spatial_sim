#include"ebola.h"
#include<ctime>
#include<limits>

using namespace std;

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
	if (argc<intNoArgs+1)
	{
		SR::srerror("First two arguments parameter_file_name and output_file_name. Rest parsed as parameter values.\n");
	}

	string strParamFile, strOutputFile, strHouseholdFile, strArgs, strHouseholdAgeDistributionFile;
	strParamFile = argv[1];
	strOutputFile = argv[2];

	if ((argc-(intNoArgs+1))%3!=0)
	{
		SR::srerror("An even number of parameter arguments are required.\n");
	}

	for (int i=intNoArgs+1;i<argc;i=i+3)
	{
		strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";
	}

	SR::ParameterSet ukPars;
	ukPars.ReadParamsFromFile(strParamFile);
	if (argc > 4)
	{
		ukPars.ReadParams(strArgs);
	}

	// Set up the population age distribution file
	strHouseholdAgeDistributionFile = ukPars.GetTag("strHouseholdAgeDistributionFile");

	// Set up the population density file
	strHouseholdFile = ukPars.GetTag("strHouseholdDensityFile");
	SR::DensityField PopulationDensityField(strHouseholdFile);
	ukPars.AddValue("dblXGridMin",PopulationDensityField.GetMinX());
	ukPars.AddValue("dblXGridSize",PopulationDensityField.GetMaxX()-PopulationDensityField.GetMinX());
	ukPars.AddValue("dblYGridMin",PopulationDensityField.GetMinY());
	ukPars.AddValue("dblYGridSize",PopulationDensityField.GetMaxY()-PopulationDensityField.GetMinY());
	SetUpExtraParameters(ukPars);

	// Initialise the ukPars.RunP.intSeed
	ukPars.intSeed = -1*ukPars.intSeed;
	gsl_rng_set(glob_rng, ukPars.intSeed);

	// Initialise Hexagon
	SR::Hexagon tmphex(ukPars);

	SR::GridHex * ukGridHexGenerate;

	//Make Grid
	ukGridHexGenerate = new SR::GridHex(ukPars,tmphex, PopulationDensityField, strHouseholdAgeDistributionFile);

	//Save Grid
	cerr << "Writing GridHex...\n";
	ofs.open((strOutputFile+"_initial_gridhex.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs << *ukGridHexGenerate;
	ofs.close();
	cerr << "...done\n";

	//Save Parameters
	cerr << "Writing Parameters...\n";
	ofs.open((strOutputFile+"_initial_params.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs << ukPars;
	ofs.close();
	cerr << "...done\n";

	gsl_rng_free (glob_rng);
	ukGridHexGenerate->~GridHex();
	return 0;

}
