#ifndef SR_INC_GRIDHEX
#define SR_INC_GRIDHEX

#include<iostream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
// #include"nr.h"
#include"SR_Utility.h"
#include"BytePackedInt.h"
#include"SR_DensityField.h"
#include"SR_PagesForThings.h"

/*

S Riley

*/

using namespace std;

namespace SR {

	class Hexagon;
	class GridHex;
	class NodeMask;

	class Node {
		friend class Hexagon;
		friend class GridHex;
	private:
		static GridHex* ptGridHex;
		float dblX,dblY;			// As is
		int ptHexInt;				// As is
		int ptSelfInt;				// As is
		PagedThingPointer<SR::Node> ptFirstHouseholdMember;
#ifdef SR_BYTEPACKED
		BytePackedInt bp1,bp2;
		void IncrementBitwiseCharacteristic() {bp1.Increment(2,7);};
		void DecrementBitwiseCharacteristic() {bp1.Decrement(2,7);};
#else
		int intVaccinationClass;	// 1 bit
		int intCharacteristic;		// 5 bits
		int generation;			// 6 bits
		int no1,no2;			// 8 bits each
		int intCurrentPointer;		// 8 bits
		int intNoLevelsQuarantine;	// 3 bits
		int intKernelIndex;		// 4 bits for now
		float fltContactAverage;
		bool blContactsFlag;
                int intAge;
		void IncrementBitwiseCharacteristic() {intCharacteristic++;};
		void DecrementBitwiseCharacteristic() {intCharacteristic--;};
#endif
	public:
#ifdef SR_BYTEPACKED
		inline void SetNoContacts(float h) {bp2.Set(18,32,static_cast<int>(h*100));};		// Contacts should only be used when current pointer isn't
		inline float GetNoContacts() {return static_cast<float>(bp2.Get(18,32))/100+0.005;};
		inline bool GetContactsFlag() {return static_cast<bool>(bp2.Get(17,17));};
		inline void SetContactsFlag() {bp2.Set(17,17,1);};
		inline void UnSetContactsFlag() {bp2.Set(17,17,0);};
		inline void SetVaccinationClass(int i) {bp1.Set(1,1,i);};
		inline int GetVaccinationClass() {return bp1.Get(1,1);};
		inline int GetCharacteristic() {return bp1.Get(2,7);};
		inline void SetCharacteristic(int c) {bp1.Set(2,7,c);};
		inline void SetGeneration(int g) {bp1.Set(8,13,g);};
		inline int GetGeneration() {return bp1.Get(8,13);};
		inline void IncrementHouseholdMax() {bp2.Increment(1,8);};
		inline void IncrementNoSpatialNeighbour(){bp2.Increment(9,16);};
		inline void SetHouseholdMax(int h) {bp2.Set(1,8,static_cast<unsigned int>(h));};
		inline void SetSpatialMax(int h) {bp2.Set(9,16,h);};
		inline void DecrementHouseholdMax() {bp2.Decrement(1,8);};
		inline void DecrementNoSpatialNeighbour() {bp2.Decrement(9,16);};
		inline int GetHouseholdMax() {return static_cast<int>(bp2.Get(1,8));};
		inline int GetNoSpatialNeighbour() {return static_cast<int>(bp2.Get(9,16));};
		inline int GetCurrentPointer() {return static_cast<int>(bp2.Get(25,32));};
		inline void IncrementCurrentPointer() {bp2.Increment(25,32);};
		inline void DecrementCurrentPointer() {bp2.Decrement(25,32);};
		inline void SetCurrentPointer(int v) {bp2.Set(25,32,static_cast<unsigned int>(v));};
		inline void IncrementQuarantineLevel() {bp1.Increment(14,18);};
		inline void DecrementQuarantineLevel() {bp1.Decrement(14,18);};
		inline void ZeroQuarantineLevel() {bp1.Set(14,18,0);};
		inline int GetQuarantineLevel() {return static_cast<int>(bp1.Get(14,18));};
		inline void SetKernelIndex(int i) {bp2.Set(19,22,i);};
		inline void GetKernelIndex(int i) {bp2.Get(19,22,i);};
#else
		inline void SetNoContacts(float h) {fltContactAverage=static_cast<float>(h-std::fmod(h,float(0.01))+0.005);};
		inline float GetNoContacts() {return fltContactAverage;};
		inline void SetAge(int a) {intAge=static_cast<int>(a);};
		inline int GetAge() {return intAge;};
		inline bool GetContactsFlag() {return blContactsFlag;};
		inline void SetContactsFlag() {blContactsFlag=1;};
		inline void UnSetContactsFlag() {blContactsFlag=0;};
		inline void SetVaccinationClass(int i) {intVaccinationClass=i;};
		inline int GetVaccinationClass() {return intVaccinationClass;};
		inline int GetCharacteristic() {return intCharacteristic;};
		inline void SetCharacteristic(int c) {intCharacteristic=c;};
		inline void SetGeneration(int g) {generation=g;};
		inline int GetGeneration() {return generation;};
		inline void IncrementHouseholdMax() {no1++;};
		inline void IncrementNoSpatialNeighbour(){no2++;};
		inline void SetHouseholdMax(int h) {no1=h;};
		inline void SetSpatialMax(int h) {no2=h;};
		inline void DecrementHouseholdMax() {no1--;};
		inline void DecrementNoSpatialNeighbour() {no2--;};
		inline int GetHouseholdMax() {return no1;};
		inline int GetNoSpatialNeighbour() {return no2;};
		inline int GetCurrentPointer() {return intCurrentPointer;};
		inline void IncrementCurrentPointer() {intCurrentPointer++;};
		inline void DecrementCurrentPointer() {intCurrentPointer--;};
		inline void SetCurrentPointer(int v) {intCurrentPointer=v;};
		inline void IncrementQuarantineLevel() {intNoLevelsQuarantine++;};
		inline void DecrementQuarantineLevel() {intNoLevelsQuarantine--;};
		inline void ZeroQuarantineLevel() {intNoLevelsQuarantine=0;};
		inline int GetQuarantineLevel() {return intNoLevelsQuarantine;};
		inline void SetKernelIndex(int i) {intKernelIndex=i;};
		inline int GetKernelIndex() {return intKernelIndex;};
#endif
		Node(int i,double x, double y);
		Node();
		double Distance(SR::Node* ptN);
		inline double GetX() {return dblX;};
		inline double GetY() {return dblY;};
		int GetIndex();
		void AddToNeighbourList(Node* ptN1);
		SR::Node** GetFirstHouseholdMember();
		inline void DecrementCharacteristic();
		void MakeCharacteristicEqualTo(int nc);
		SR::Node* GetMaximalNode();
		Hexagon* GetHexagon();
		string OutputIndexLocationNetwork();
		string OutputNodeAndHexagonDetails();
		inline void SetX(double x) {dblX=x;};
		inline void SetY(double y) {dblY=y;};
		friend ofstream& operator<<(ofstream& ofs, Node& n);
		friend SR::Node ReadNodeBinaryFromFile(ifstream& ifs);
		string OutputNodeToLine();
		inline GridHex* GetPtGridHex(){return ptGridHex;};
		inline void SetPtGridHex(GridHex* g) {ptGridHex = g;};
		inline void SetPtSelf(SR::Node** pt);
		inline SR::Node** GetPtSelf();
	};

	// Hexagon (not abstract)
	class Hexagon {
		friend class Node;
		friend class GridHex;
	private:
		static const int maxcb = 15;
		static const int maxmn = 2;
		bool blMembersAltered;
		double dblCenX, dblCenY;
		int intCoordX, intCoordY, usedcb, usedmn;
		int ptptFirstNode;
		int ptptLastNode;
		int intRegionalTreatmentStatus;
		GridHex *ptGrid;
		int arrIntCharBoundaries[maxcb];
		SR::Node arrMaximalNodes[maxmn];
		void AddNode(SR::Node* ptNode);
		void IncrementCharOfNode(SR::Node** ptptNode, int schar);
		void DecrementCharOfNode(SR::Node** ptptNode, int schar);
		void AssignCentrePoint(double hw, SR::ParameterSet& p);
	public:
		inline double GetX() {return dblCenX;}
		inline double GetY() {return dblCenY;}
		inline int GetCoordX() {return intCoordX;}
		inline int GetCoordY() {return intCoordY;}
		void SetCoordX(int x) {intCoordX=x;}
		void SetCoordY(int y) {intCoordY=y;}
		Hexagon(){};
		Hexagon(SR::ParameterSet& p);
		Hexagon(int chs, int mn);
		SR::Node** FirstOfChar(int i);
		SR::Node** LastOfChar(int i);
		string OutputAllIndexAndLocations();
		string OutputIndexByCharacteristicOnOneLine();
		void SetMaximalNode(SR::Node* vn, int index);
		SR::Node* GetMaximalNode(int index);
		double Distance(SR::Node* ptN);
		SR::Node** GetFirstNode();
		SR::Node** GetLastNode();
		int GetNoNodes();
		void SetFirstNode(SR::Node** ptpt);
		void SetLastNode(SR::Node** ptpt);
		inline GridHex* GetPtGrid() {return ptGrid;};
		inline void SetPtGrid(GridHex* g) {ptGrid=g;};
		friend ofstream& operator<<(ofstream& ofs, Hexagon& h);
		friend Hexagon ReadHexagonBinaryFromFile(ifstream& ifs);
		string OutputHexagonToMultipleLines();
		inline void IncrementLastNode() {ptptLastNode++;};
		inline void IncrementFirstNode() {ptptFirstNode++;};
		inline int GetRegionalTreatment() {return intRegionalTreatmentStatus;};
		inline void SetRegionalTreatment(int s) {intRegionalTreatmentStatus = s;};
		inline void SetMembersAltered() {blMembersAltered=true;};
		inline void UnSetMembersAltered() {blMembersAltered=false;};
		inline bool GetMembersAltered() {return blMembersAltered;};
		void ResetMembersCharGen(int c, int g);
		int GetNoInfectedNotInfectious();
	};

	// Grid of hexagons
	class GridHex {
	private:
		double dblHexagonWidth;
		int intMinXCoord,intMinYCoord,intNoXCoords,intNoYCoords;
		double maxdx, maxdy;
		SR::Hexagon *lastHexagon;
		PagesForThings<SR::Node>* Blocks;
		SR::Node *vecNodes;
		int sizeVecNodes;
		SR::Hexagon *vecHexagon;
		int sizeVecHexagons;
		SR::Node **vecPtNodesHexOrder;
		IntCoord RealToHexCoords(double x, double y);
	public:
		GridHex(SR::ParameterSet& p, Hexagon tmphex, DensityField& houses);
		GridHex(SR::ParameterSet& p, Hexagon tmphex, ifstream& ifs);
		~GridHex();
		void ReserveMemoryForNetwork(SR::PagesForThings<SR::Node>* blks);
		int CalculateSizeOfNetwork();
		bool WriteNodeLocationsToFile(string filename);
		bool WriteNodeLocationsAndSizesToFile(string filename);
		bool WriteArcsToFile(string filename);
		bool WriteNodeLocationsAndLinksToFile(string filename);
		inline SR::Node* FirstNode() {return vecNodes;};
		inline SR::Node* LastNode() {return vecNodes+sizeVecNodes;};
		inline SR::Hexagon* FirstHexagon() {return vecHexagon;};
		inline SR::Hexagon* LastHexagon() {return lastHexagon;};
		inline int GetNoHexagons() {return LastHexagon()-FirstHexagon();};
		void HouseholdSetupRoutine(vector<vector<int> >& vn, UNTIMEDEVENT ue);
		void MakeAllNodesThisCharacteristicAndGeneration(int c, int g);
		string OutputCharacteristicsOnOneLine();
		string OutputHexagonCharacteristicsOnManyLines();
		double GetHexagonWidth() {return dblHexagonWidth * 6378.7;}; // Max width of a hexagon a hex (at the equator)
		double CalculateAverageSpatialNeighbours();
		double CalculateAverageHousehold();
		inline int GetNoNodes() {return sizeVecNodes;};
		void AssignHouseholds();
		inline PagesForThings<SR::Node>* GetPtBlocks() {return Blocks;};
		inline SR::Node** GetFirstOfVecPtsHexagon() {return vecPtNodesHexOrder;};
		friend ofstream& operator<<(ofstream& ofs, GridHex& gh);
		string OutputGridHexToMultipleLines();
		inline double GetMaxDx() {return maxdx;};
		inline double GetMaxDy() {return maxdy;};
		void AssignToPagesOfNodes(PagesForThings<SR::Node>* blks);
		double CalcExpectedSpatial(SR::ParameterSet& p, KERNEL k);
		void MakeAllHexagonsNotInRegionalTreatment();
		int GetTotalInfectedNotInfectious();
	};

	class NodeMask {
	private:
		bool *mask;
		int  *vecNodesSeen;
		int maxNoNodes, seenNoNodes;
	public:
	    NodeMask(){SR::srerror("No default constructor for NodeMask");};
	    NodeMask(SR::GridHex &GH);
	    ~NodeMask();
	    void AgeMask(int lbInc, int ubInc, SR::GridHex &GH);
	    void NullMask(SR::GridHex &GH);
	    inline int GetNoSeenNodes(){return seenNoNodes;};
	};

	void AssignHexagons(SR::Hexagon* ptHexStart, SR::Hexagon* ptHexEnd, GridHex &g);
	Node ReadNodeBinaryFromFile(ifstream& ifs);
	Hexagon ReadHexagonBinaryFromFile(ifstream& ifs);

};

#endif
