#ifndef SR_INC_WORKPLACES
#define SR_INC_WORKPLACES

#include<iostream>
// #include"nr.h"
#include"SR_Utility.h"
#include"SR_GridHex.h"
#include"SR_EventMatrix.h"
#include"SR_Kernels.h"
#include"CachedLookups.h"


/*

[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {

	struct workplace {
		double x, y;
		int index, firstnode, totalnodes, currentlyusednodes;
	};

	class Workplaces {
	private:
		double startscale; // constant! initialised in contructor
		static const int noOrders = 7;
		static const int incsPerOrder = 50;
		static const int outputBinNumber = 2000;
		int* vecNodeIntsInWorkplaceOrder;
		int sizeVecNodeIntsInWorkplaceOrder;
		int* vecWorkplaceIntsForEachNode;
		int sizeVecWorkplaceIntsForEachNode;
		workplace* vecWorkplaces;
		int sizeVecWorkplaces;
		int* vecDistFrequencies;
		int sizeVecDistFrequencies;
		int* vecDistFrequenciesOld;
		double outputBinSize;
		int WorkplaceSize(int i);
		void MoveNode(int nd, int wkf);
		void MoveNodeDash(int nd, int wkf);
		void MakeWorkplacesConsistentWithNodes(SR::GridHex& g, SR::ParameterSet& p);
		inline int GetDistanceBin(double d);
		int intNoNodes,intNoWorkplaces,intNoWorkPlaceGrids,intNoXWorkplaceGrids,intNoYWorkplaceGrids;
		workplace **ptFirstOrderedWorkplace;
		int *ptIntFirstNoWorkplacesInGrid;
		workplace ***ptFirstPtWorkPlace;
		int *ptArrNumberNearbyWorkplaces;
		workplace ***ptArrPtPtNearbyWorkplaces;
		int *ptArrNumberEachNearbyWorkplaces;
		double dblGridSize;
		inline pair<int,int> CalcGridXY(double x, double y);
		inline bool GridIsClose(pair<int,int>& a, pair<int,int>& b);
		inline int GetGridIndex(pair<int,int>& p);
		void ArchiveCommuteDist();
		int degreesfreedom;
	public:
		Workplaces(SR::GridHex& g, SR::ParameterSet& p, SR::DensityField& wpd);
		Workplaces() {SR::srerror("No default constructor for SR::Workplaces");};
		~Workplaces();
		void WriteWorkplacesToFile(string f);
		void WriteCommutesToFile(SR::GridHex& g, string f);
		inline workplace* GetWorkplaceOfNode(int i);
		void MCMCUpdate(GridHex& gh, ParameterSet& p,double pdistance(double,int,double*), string funcfile, int minSamplesMillions, int maxSamplesMillions, int& cn);
		void GenerateDistributionOfPossibleCommutes(GridHex& gh, ParameterSet& p, double grid_dx, double hist_dx, string funcfile);
		int* FirstCollegue(int n);
		int* OnePastLastCollegue(int n);
		int* VeryFirstCollegue() {return vecNodeIntsInWorkplaceOrder;};
		int* OnePastVeryLastCollegue() {return vecNodeIntsInWorkplaceOrder+sizeVecNodeIntsInWorkplaceOrder;};
		string PrintDistributionOfCommutes();
		string PrintDistributionOfSizes();
		pair<double,double> CalcChiSquared(double zed_val, int df);
	};

	void GenerateAllNeighbours(GridHex& gh, ParameterSet& p, PROCESS proc1, KERNEL kern1, KERNEL kern2, UNTIMEDEVENT ue, Workplaces &w);

}

#endif
