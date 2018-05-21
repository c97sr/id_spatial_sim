#ifndef SR_INC_UK2K
#define SR_INC_UK2K

// STL Includes
#include<iostream>
#include<string>
#include<bitset>

// Third party includes
#include"nr.h"

// My inlcudes
#include"SR_CachedLookups.h"
#include"SR_GridHex.h"
#include"SR_Workplaces.h"
#include"SR_InitialConditions.h"
#include"SR_GslRng.h"

using namespace std;

bool evCountHouseholdNeighbour(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evCountSpatialNeighbour(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evAddToNeighbourList(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evInfection(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evBecomeInfectious(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterEarlyRash(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterLateRash(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evRecover(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evDie(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evRecover(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evVaccinate(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evBecomeLatentVaccinatedTwo(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterContactQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterFeverQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterRashQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEndQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
bool evEnterLateRash(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
void procSpatialNeighbourSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procSpatialNeighbourAdd(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procSubsequentToInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procContactTrace(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procRegionalVaccination(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procEnterProdrome(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
void procEnterQuarantine(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
double kernNeighbourSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernSpatialInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernSpatialInfectionOff(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernFileCached(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernFileCachedEarlyRash(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernHouseholdInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernNeighbourInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernHouseholdInfectionSymp(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernNeighbourInfectionSymp(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernIntroSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double kernSpatialEbola(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
double DensityTest(double x, double y);
void SetUpExtraParameters(SR::ParameterSet &p);
void RunExtraParameters(SR::ParameterSet &p);
void CalcBetasWithAttackRate(SR::ParameterSet &p, SR::GridHex& g, SR::KERNEL k);
double CalcProbContactInfected(double mean, int alpha, double timestep, double intstep, int maxstep, double hazard);
vector<double> GenerateGammaVector(double mean, int alpha, double timestep, double intstep, int maxstep);
double CalcProbContactInfected(vector<double>& gamvec, double timestep, double hazard);
double fnLogLogBiphasic(double d, SR::ParameterSet& p);
double fnSquaredLogLogBiphasic(double d, SR::ParameterSet& p);
double fnOffsetPower(double d, SR::ParameterSet& p);
double fnOffsetPower(double d, int np, double* p);
double distTestDist(double ave, SR::ParameterSet &p);
int PoxEventToInt(SR::UNTIMEDEVENT ue);
SR::UNTIMEDEVENT PoxIntToEvent(int e);
void setGZRegionMembership(SR::GridHex &g);

class GatherUKPoxInfoA : public SR::GatherRunInformation {
	struct logEvent {
		double time,x,y,infectx,infecty;
		int nodeindex,realisation,eventindex,generation,infectorindex;
	};
	struct coords {
		double x,y;
	};
private:
	static const int maxstringlength = 200;
	static const int noEvents=8;
	int noTimeSteps,noRealisations;
	int noRealisationsCompleted;
	int intCurrentTimestep;
	int intSizeFileBase;
	double maxSpatialExtent;
	double dblTimeStep;
	double cumAverage;
	int intCurrentStackSize;
	int intCurrentSpatialStackSize;
	int intLastSpatialStackSize;
	logEvent tmpLogEvent;
	vector<int> vecInc;
	vector<double> processData;
	vector<double> vecTimes;
	vector<double> vecAverages;
	vector<double> vecStandardDeviations;
	vector<double> vecCumTotals;
	vector<logEvent> eventStack;
	vector<coords> locationStack;
	double cumStandardDeviation;
	char fileBase1[maxstringlength];
	void SetFileBase(string s);
	bool blSpatialMeasuresOn;
	bool blLogEvents;
public:
	string GetFileBase();
	bool FlushStack();
	void ResetRealisations() {noRealisationsCompleted=0;};
	void ResetTimestep() {intCurrentTimestep=0;maxSpatialExtent=0;intCurrentSpatialStackSize=0;intLastSpatialStackSize=0;};
	inline void IncrementRealisations() {noRealisationsCompleted++;};
	void IncrementTimestep();
	vector<double>::iterator AccessData(int ne, int nt, int nr);
	GatherUKPoxInfoA(int nt, double dt, int nr, string f, int ss, int lss);
	GatherUKPoxInfoA(){SR::srerror("No default constructor for GatherUKPoxInfoA");};
	void RegisterEventAfterApplication(vector<SR::UntimedEvent>::iterator ptEv, int cb);
	string ConvertEventPointerToString(SR::UNTIMEDEVENT ue);
	string CurrentCumPrev();
	string OutputRealisationData(int index);
	void ClearCumPrev();
	void UpdateIncidenceFile(int index);
	void CalcForIndex(int index);
	string OutputTable(int index);
	string OutputTableForIDL(int index);
	bool WriteAllEventIncidencesToFile();
	void UpdateAtEndOFTimeStep();
	bool AddNodeLocationToSpatialStack(SR::Node *ptN);
	void RecalculateSpatialExtent();
	bool AreSpatialMeasuresOn(){return blSpatialMeasuresOn;};
	void SpatialMeasuresOff() {blSpatialMeasuresOn=false;};
	bool EventStackFull();
	void StartEventLogging() {blLogEvents=true;};
	void StopEventLogging() {blLogEvents=false;};
};

#endif
