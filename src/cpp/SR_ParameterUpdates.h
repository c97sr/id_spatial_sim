#ifndef SR_INC_PARAMETERUPDATES
#define SR_INC_PARAMETERUPDATES

#include"SR_Parameter.h"
#include"SR_Utility.h"

using namespace std;

namespace SR {
	class ParameterUpdates {
	private:
		static const int mp = 30;
		Label outputFile,inputFile,uniqueId;
		int noParams,sampleNo,intNoChosenParam;
		double currentLnLike,dblStoredParam,dblLogProposalFactor;
		double totalWeight;
		ParameterSet* ps;
		// int noParamsFitted;
		vector<Label> paramNames;
		vector<double*> paramPointers;
		vector<double> paramInitial;
		vector<double> paramCurrent;
		vector<double> paramWeight;
		vector<double> paramMin;
		vector<double> paramMax;
		vector<double> paramPropStep;
		vector<bool> paramRangeIsLog;
	public:
		ParameterUpdates() {SR::srerror("No default constructor.");};
		ParameterUpdates(SR::ParameterSet& ps_in, string inFileName, string outFileName, string uniqueid);
		int GetNoParameters() {return noParams;};
		string GetParameterName(int n) {return paramNames[n].Get();};
		double GetParamValueByIndex(int i) {return *(paramPointers[i]);};
		double GetProposalFactor() {return dblLogProposalFactor;};
		int GetNoSamples() {return sampleNo;};
		void ProposeUpdate();
		void AcceptUpdate();
		void DeclineUpdate();
		void LogParameterValues(double lnlike, int count);
		void LogParameterValues(double lnlike, ostringstream &oss, int count);
	};
}

#endif
