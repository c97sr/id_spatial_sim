#ifndef SR_INC_KERNELS
#define SR_INC_KERNELS

#include<iostream>
// #include"nr.h"
#include"SR_Utility.h"
#include"SR_GridHex.h"
#include"SR_Stats.h"
#include"SR_EventMatrix.h"

/*

Defined the partially abstract kernel class and then the various specific kernel classes derived from it.

sr@stevenriley.net

*/

using namespace std;

namespace SR {

	class Kernel {
	protected:
		vector<int> vecIntSourceCharacteristics;
		vector<int> vecIntTargetCharacteristics;
		PROCESS procProcess;
		KERNEL kernKernel;
		int intMaximalIndex;
		double dblExpNumber;
		double dblCurrentConst;
		virtual void GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode);
		void GenerateEventsOneSourceNodeWithCross(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode, SR::ONE_D_DISTRIBUTION c, double* ptC);
	public:
		Kernel(vector<int>& sc, vector<int>& tc, KERNEL k, PROCESS p, int mi);
		void GenerateAllEvents(SR::ParameterSet& p, EventMatrix& em, GridHex& gh);
		void GenerateAllEventsWithCross(SR::ParameterSet& p, EventMatrix& em, GridHex& gh, SR::ONE_D_DISTRIBUTION c, double* ptC);
		virtual void OneNodeOneHex(EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double cp);
	};

	// Hex kernel from here

	class HexKernel { // note hex kernel not based on Kernel above
	protected:
		vector<int> vecIntSourceCharacteristics;
		SR::Hexagon **vecPtHex;
		SR::Hexagon **ptOnePastLast;
		PROCESS procProcess;
		double dblMaxVax;
		int	intPopSize;
		int globalMaxNoTreatments;
		int maxHex;
		double globalMaxTreatmentRate;
		double localMaximumRate;
		void GenerateEventsForOneHex(SR::Hexagon *ptHex, SR::EventMatrix &em, SR::ParameterSet &p);
	public:
		HexKernel() {SR::srerror("No default constructor for HexKernel");};
		~HexKernel();
		HexKernel(vector<int>& sc, PROCESS proc, int mh, double gmr, int gmn, double lmr);
		void GenerateEvents(EventMatrix &em, ParameterSet& p);
		void AddHexagons(GridHex &gh, double x, double y, double range);
		void AddHexagons(GridHex &gh, SR::Hexagon *ptHex, double range);
		void RemoveHexagon(SR::Hexagon *ptHex);
		void UpdateRegion(GridHex &gh, double range, double max_symp, double max_sus);
		void ResetRegion(GridHex &gh);
	};

	// Spatial kernel from here
	class SpatialKernel : public Kernel {
	private:
		int *ptStack, *ptStackTop;
		int i,intTotalPoss, intNumberChosen, intNodeToTry;
		bool boolLoopTest;
		double probability;
		int* vecStackOfInts;
		int *ptInt, *ptIntEnd;
		int sizeofstack;
		int intLevelOfStack;
		void AddToStack(int i);
		void GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode);
	public:
		SpatialKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int ss, int mi);
		~SpatialKernel();
		void EagerOneNodeOneHex(EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p);
		void LazyOneNodeOneHex(EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double cp);
		void OneNodeOneHex(EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double cp);
	};

	// Household kernel from here

	class HouseholdKernel : public Kernel {
	protected:
		void GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode);
	public:
		HouseholdKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int mi);
	};

	// Neighbour kernel from here

	class NeighbourKernel : public Kernel {
	protected:
		void GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode);
	public:
		NeighbourKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int mi);
	};

}

#endif
