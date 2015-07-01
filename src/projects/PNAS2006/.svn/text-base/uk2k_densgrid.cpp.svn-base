#include"uk2k.h"

int main(int argc, char* argv[]) {
#ifdef _DEBUG 
	cerr << "DEBUG version" << endl; 
#endif
	if (argc!=9) SR::srerror("Density file, output file, minx, maxx, gapx, miny, maxy, gapy required as arguments.\n");
	string strArgs, strOutputFile, strDensityFile;
	double xmin, xmax, xgap, ymin, ymax, ygap;
	for (int i=1;i<argc;++i){strArgs+=argv[i];strArgs+="\t";};
	istringstream iss(strArgs);
	iss >> strDensityFile >> strOutputFile >> xmin >> xmax >> xgap >> ymin >> ymax >> ygap;
	cerr << "Generating density field from " << strDensityFile << "\n";
	SR::DensityField df(strDensityFile, xgap, ygap, xmin, xmax, ymin, ymax);
	ofstream ofs(strOutputFile.c_str());
	if (ofs.fail()) SR::srerror("Problem with output file");
		ofs << df.Table();
	ofs.close();
	return 0;
}