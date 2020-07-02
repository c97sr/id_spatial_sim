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


	// Read from command line and set up parameter object
	int intNoArgs = 2;
	if (argc<intNoArgs+1)
	{
		SR::srerror("First two arguments parameter_file_name and output_file_name.\n");
	}

	string strParamFile, strOutputFile, strInputFile;
	strParamFile = argv[1];
	strInputFile = argv[2];
	strOutputFile = argv[3];

	//Load Grid
	cerr << "Loading GridHex...\n";
	ifs.open((strInputFile+"_initial_gridhex.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary gridhex file.");
	SR::GridHex * ukGridHex;
	ukGridHex = new SR::GridHex();
	ifs >> *ukGridHex;
	ifs.close();
	cerr << "...done\n";
	cerr << ukGridHex->GetNoNodes();

	//Load Parameters
	ifs.open((strInputFile+"_initial_params.hex").c_str(),ios::binary);
	if (ifs.fail()) SR::srerror("Problem opening binary parameter file.");
	SR::ParameterSet ukPars;
	ifs >> ukPars;
	ifs.close();

	ukPars.ReadParamsFromFile(strParamFile);

	//Generate Masks
	cerr << "Generating masks...\n";

	string strMaskProportionsFile = ukPars.GetTag("strStrataProportionsFile");
	ifs.open(strMaskProportionsFile);
	if (ifs.fail()) SR::srerror("You idiot.");
	SR::GroupDistribution* dist = new SR::GroupDistribution(ifs, ukGridHex);
	ifs.close();
	cerr << "...done\n";

	//Save Masks
	for(int i = 0; i < dist->GetNoMasks(); i++) {
		cerr << "Saving mask " << i + 1 << "\n";
		ofs.open((strOutputFile+"_group_distribution_"+to_string(i)+".hex").c_str(),ios::binary);
		if (ofs.fail()) SR::srerror("You idiot.");
		dist->WriteMask(ofs, i);
		ofs.close();
		cerr << "done\n";
	}

	//Save Parameters

	ofs.open((strOutputFile+"_initial_params.hex").c_str(),ios::binary);
	if (ofs.fail()) SR::srerror("You idiot.");
	ofs << ukPars;
	ofs.close();

	ukGridHex->~GridHex();
	dist->~GroupDistribution();

	return 0;
}

