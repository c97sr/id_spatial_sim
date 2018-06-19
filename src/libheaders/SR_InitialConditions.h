#ifndef SR_INC_INITIALCONDITIONS
#define SR_INC_INITIALCONDITIONS

#include<iostream>
// #include"nr.h"
#include"SR_Utility.h"
#include"SR_GridHex.h"

/*
[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {

	class InitialConditions {
	private:
		int* ptvecListOfHexagons;
		int* ptvecNodesInHexagons;
		int intTotalNodes;
		int intHexagonsUsed;
		vector<int> vecNodes;
		double x,y,r;
		int maxsamples,currentsamples,initialnumber;
		bool blLocalSeeding;
	public:
		InitialConditions() {SR::srerror("No default constructor for InitialConditions.");};
		InitialConditions(GridHex& g, ParameterSet& p, int no, int msn, double x_in, double y_in, double r_in, int ms);
		~InitialConditions();
		void Reselect(GridHex& g, ParameterSet& p);
		void Reselect(GridHex& g, ParameterSet& p, int number);
		void ApplySeed(GridHex& g, SR::ParameterSet& p, PROCESS proc, EventMatrix& em);
		void Trickle(GridHex& g, SR::ParameterSet& p, PROCESS proc, EventMatrix& em, double trate, double tdur);
	};

}

#endif
