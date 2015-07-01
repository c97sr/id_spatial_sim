#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include<algorithm>
#include<iomanip>
#include<time.h>

using namespace std;

namespace SR {

	// Rate parameter structures
	class RateParam {
	private:
		static const int maxstringlength = 100;
		double dblVal,dblMin,dblMax;
		bool blFitted,blLog;
		char strName[maxstringlength];
		int nameSize;
	public:
		RateParam() {};
		RateParam(ifstream &ifs);
		void SetStringAndInitialValue(string s,double val);
		inline void ChangeValue(double nv) {dblVal = nv;}; // This can have bound checking eventually
		inline void ChangeMax(double nv) {dblMin = nv;};
		inline void ChangeMin(double nv) {dblMax = nv;};
		inline double GetValue() {return dblVal;};
		inline double GetMax() {return dblMin;};
		inline double GetMin() {return dblMax;};
		inline double* GetPtValue() {return &dblVal;};
		inline double* GetPtMax() {return &dblMax;};
		inline double* GetPtMin() {return &dblMin;};
		inline string GetName() {string rtnstr(strName,nameSize); return rtnstr;};
		friend ofstream& operator<<(ofstream& ofs, RateParam& rp);
		friend RateParam ReadRateParamBinaryFromFile(ifstream& ifs);
	};

	// Run parameter structures
	struct RunParam {
		int intNoRuns,intNoParamSamples;
		int intSeed;
		int intNoNodes;
		int intNoHousehold;
		int intMaxDailyEvents,intMaxDelayTimesteps;
		int intRealisationsPerParameterSet,intCurrentRealisation,intInfectionKernelStackSize;
		int intNoCharacteristics, intNoMaximalNodes;
		int intBlockSize, intNumberOfBlocks;
		int intOutputPrecision,intIndexOfSusceptible;
		int intMaxNeighbourEvents, intIndexOfMaximalElement;
		int intMaxNoHexagons;
		double dblTimeStep,dblStartTime,dblEndTime,dblCurrentTime;
		double dblXGridSize,dblYGridSize;
		double dblMaxSpatialNeighbours,dbldblMaxWithinHexSampleRate,dblHexagonWidth;
		double dblAverageDensity,dblAverageHousehold;
		int intNumberOfWorkplaces;
		double dblPropColleguesInNetwork;
		int intNoMcmcSteps;
		double dblXGridMin, dblYGridMin;
		double dblAverageWorkplaceSize;
	};

		class Parameters {
	private:
		bool boolLockedForPointers;
		int intNumberRateParametersUsed;
		vector<RateParam>::iterator FindRateParam(string s);
		vector<RateParam>::iterator FindRateParamNoError(string s);
	public:
		RunParam RunP;
		vector<RateParam> RateP;
		Parameters() : RateP(0) {SR::srerror("Cannot call default constructor for parameters.");}; 
		Parameters(int n) : RateP(n) {boolLockedForPointers=false;intNumberRateParametersUsed=0;}; 
		Parameters(ifstream &ifs, int n); 
		void ReadFromFile(string f, string fname);
		void AddRateParameter(string s, double val);
		inline void Lock() {boolLockedForPointers=true;}; // no "unlock" defined!
		void ChangeRateParameterValue(string s, double d);
		double ParamValue(string s);
		double* ParamPointer(string s);
		void AssignRateParamsFromString(string s);
		string OutputAllRateParamValues();
		friend ofstream& operator<<(ofstream& ofs, Parameters& p);
		friend Parameters ReadParametersBinaryFromFile(ifstream& ifs);
		bool RateParamDefined(string s);
	};

}

void SR::Parameters::AddRateParameter(string s, double val) {	
	if (boolLockedForPointers) SR::srerror("Parameters cannot be added when parameter structure locked.\n");
	if (FindRateParamNoError(s)!=RateP.begin()+intNumberRateParametersUsed && intNumberRateParametersUsed>0) SR::srerror("Attempting to add duplicate parameter name in SR::Parameters::AddRateParameter.\n");
	if (intNumberRateParametersUsed==RateP.size())SR::srerror("Maximum number of rate parameters already reached");
	SR::RateParam tmpRp;
	tmpRp.SetStringAndInitialValue(s,val);
	RateP[intNumberRateParametersUsed]=tmpRp;
	intNumberRateParametersUsed++;
};

void SR::Parameters::ReadFromFile(string s, string fname) {
	if (boolLockedForPointers) SR::srerror("Parameters cannot be added when parameter structure locked.\n");
	if (FindRateParamNoError(s)!=RateP.begin()+intNumberRateParametersUsed && intNumberRateParametersUsed>0) return;
	// NMF - Oh YES - fscanf rules!!!
	FILE * dat;
	char buf[255];
	double val;
	int found_flag=false;
	if(!(dat=fopen(fname.c_str(),"r"))) SR::srerror("Can't open parameter file\n");
	while((!feof(dat))&&(!found_flag)) {
		fscanf(dat,"%s",buf);
		if(strcmp(buf,s.c_str())==0){
			fscanf(dat,"%lg",&val);
			found_flag=true;
		}
	}
	fclose(dat);
	if(!found_flag) {
		cerr << "Can't find parameter '"<< s <<"' in parameter file '"<<fname<<"'\n";
		exit(1);
	}
	SR::RateParam tmpRp;
	tmpRp.SetStringAndInitialValue(s,val);
	RateP[intNumberRateParametersUsed] = tmpRp;
	intNumberRateParametersUsed++;
};

void SR::RateParam::SetStringAndInitialValue(string s,double val) {
	// if (val < dblMin || val > dblMax) SR::srerror("Parameter out of min max range when set in SetStringAndInitialValue.\n");
	if (s.size()>maxstringlength) SR::srerror("Parameter name too long in SR::RateParam::SetStringAndInitialValue");
	nameSize = s.size();
	for (int i=0;i<s.size();++i) strName[i]=s[i];
	dblVal = val;
};

vector<SR::RateParam>::iterator SR::Parameters::FindRateParam(string s) {
	vector<RateParam>::iterator rtnIt = FindRateParamNoError(s);
	if (rtnIt==RateP.begin()+intNumberRateParametersUsed) {
		cerr << "Name '" << s << "' not found in SR::Parameters::FindRateParam";
		exit(1);
	}
	return rtnIt;
};

vector<SR::RateParam>::iterator SR::Parameters::FindRateParamNoError(string s) {
	vector<RateParam>::iterator rtnIt = RateP.begin();
	vector<RateParam>::iterator locationIt = RateP.begin()+intNumberRateParametersUsed;
	while (rtnIt != RateP.begin()+intNumberRateParametersUsed) {
		if (rtnIt->GetName()==s) locationIt = rtnIt; 
		rtnIt++;
	}
	return locationIt;
};

void SR::Parameters::ChangeRateParameterValue(string s, double d) {
	vector<SR::RateParam>::iterator ptRp = FindRateParam(s);
	ptRp->ChangeValue(d);
};

double SR::Parameters::ParamValue(string s) {
	vector<SR::RateParam>::iterator ptRp = FindRateParam(s);
	return ptRp->GetValue();	
};

double* SR::Parameters::ParamPointer(string s) {
	if (!boolLockedForPointers) SR::srerror("Parameters must be locked for pointers to be assigned in SR::Parameters::ParamPointer/\n");
	vector<SR::RateParam>::iterator ptRp = FindRateParam(s);
	return ptRp->GetPtValue();	
};

void SR::Parameters::AssignRateParamsFromString(string s) {
	string name;
	double value;
	istringstream iss(s);
	while (iss>>name) {
		iss>>value;
		AddRateParameter(name,value);
	}
};

string SR::Parameters::OutputAllRateParamValues() {
	ostringstream oss;
	for (int i=0;i<intNumberRateParametersUsed;++i) {
		oss << RateP[i].GetName() << "\t" << RateP[i].GetValue() << "\n";
	}
	return oss.str();
};

SR::RateParam::RateParam(ifstream& ifs) {
	static char *filePointer;
	filePointer = (char*)(this);
	ifs.read(filePointer,sizeof(SR::RateParam));
};

ofstream& SR::operator<<(ofstream& ofs, SR::RateParam& rp) {
	static char *filePointer;
	filePointer = (char*)(&rp); ofs.write(filePointer,sizeof(SR::RateParam));
	return ofs;
};

SR::RateParam SR::ReadRateParamBinaryFromFile(ifstream& ifs) {
	SR::RateParam rtnrp;
	static char *filePointer;
	filePointer = (char*)(&rtnrp);
	ifs.read(filePointer,sizeof(SR::RateParam));
	return rtnrp;
};

ofstream& SR::operator<<(ofstream& ofs, SR::Parameters& p) {
	static char *filePointer;
	static int size;
	size = p.intNumberRateParametersUsed;
	filePointer = (char*)(&size); ofs.write(filePointer,sizeof(int));
	filePointer = (char*)(&p.boolLockedForPointers);  ofs.write(filePointer,sizeof(bool));
	filePointer = (char*)(&p.RunP); ofs.write(filePointer,sizeof(p.RunP));
	for (int i=0;i<size;++i) ofs << p.RateP[i];
	return ofs;
};

SR::Parameters SR::ReadParametersBinaryFromFile(ifstream& ifs) {
	static char *filePointer;
	static int size;
	filePointer = (char*)(&size); ifs.read(filePointer,sizeof(int));
	SR::Parameters p(size);
	filePointer = (char*)(&p.boolLockedForPointers); ifs.read(filePointer,sizeof(bool));
	filePointer = (char*)(&p.RunP); ifs.read(filePointer,sizeof(p.RunP));
	for (int i=0;i<size;++i) p.RateP[i] = SR::ReadRateParamBinaryFromFile(ifs);
	return p;
};

SR::Parameters::Parameters(ifstream &ifs, int n) : RateP(n) {
	static char *filePointer;
	static int size;
	filePointer = (char*)(&size); ifs.read(filePointer,sizeof(int));
	filePointer = (char*)(&boolLockedForPointers); ifs.read(filePointer,sizeof(bool));
	filePointer = (char*)(&RunP); ifs.read(filePointer,sizeof(RunP));
	for (int i=0;i<size;++i) RateP[i] = SR::RateParam(ifs);
};
