#include"SR_GridHex.h"

gsl_rng * glob_rng;

SR::GridHex* SR::Node::ptGridHex=NULL;

string SR::Node::OutputIndexLocationNetwork() {
	// up to here - next thing to do is to output the network of the node here... variable width lines
	ostringstream oss;
	oss << GetIndex() << "\t" << GetX() << "\t" << GetY();
	SR::Node** ptptStart = GetFirstHouseholdMember();
	oss << "\th\t";
	for (int i=0;i<GetHouseholdMax();++i) {oss << (*ptptStart)->GetIndex() << "\t";ptptStart++;}
	oss << "\ts\t";
	for (int i=0;i<GetNoSpatialNeighbour();++i) {oss << (*ptptStart)->GetIndex() << "\t";ptptStart++;}
	return oss.str();
};

void SR::Hexagon::AssignCentrePoint(double hw, SR::ParameterSet& p) {
	int x = intCoordX;
	int y = intCoordY;
	static double* xmin = p.GetPointer("dblXGridMin");
	static double* ymin = p.GetPointer("dblYGridMin");
	// Set the centre of each hexagon
	if (x < 0 || y < 0)
		SR::srerror("Grid coords not assigned.\n");
	double h=hw/2;
	double r=0.866025403784439*hw;
	if (x%2==0) {
		dblCenY = static_cast<double>(y)*2*r;
		dblCenX = static_cast<double>(x)*(hw+h);
	} else {
		dblCenY = static_cast<double>(y)*2*r+r;
		dblCenX = static_cast<double>(x)*(hw+h);
	}
	dblCenX+=*xmin;
	dblCenY+=*ymin;
};

SR::Hexagon::Hexagon(SR::ParameterSet& p) { // :
	usedcb = p.GetIntValue("intNoCharacteristics");
	usedmn = p.GetIntValue("intNoMaximalNodes");
	if (usedcb > maxcb) SR::srerror("Too many characteristics used for basic hexagon structure");
	if (usedmn > maxmn) SR::srerror("Too many maximal nodes used for basic hexagon structure");
	for (int i=0;i<usedcb;++i) arrIntCharBoundaries[i]=0;
	intCoordX=-1;
	intCoordY=-1;
	dblCenY =-1;
	dblCenX =-1;
	ptptFirstNode=0;
	ptptLastNode=0;
	SetMembersAltered();
};

SR::Hexagon::Hexagon(int chs, int mn) { // :
	usedcb = chs;
	usedmn = mn;
	if (usedcb > maxcb) SR::srerror("Too many characteristics used for basic hexagon structure");
	if (usedmn > maxmn) SR::srerror("Too many maximal nodes used for basic hexagon structure");
	for (int i=0;i<usedcb;++i) arrIntCharBoundaries[i]=0;
	intCoordX=-1;
	intCoordY=-1;
	dblCenY =-1;
	dblCenX =-1;
	ptptFirstNode=0;
	ptptLastNode=0;
	SetMembersAltered();
};

void SR::Hexagon::IncrementCharOfNode(SR::Node** ptptNode, int schar) {
	// Needs to swap the self pointers of the two nodes involved - should be a void function
	// schar is indexed from zero.
	// arrIntCharBoundaries refers to the first element of the _next_ characteristic
	//   - this approach makes best use of .end() beign one past last while .begin() is eq first
	SR::Node** tmpPtPt;
	// SR::Node* tmpPt;
	SR::Node** ptptBoundary = GetFirstNode()+arrIntCharBoundaries[schar]-1;
	swap(*ptptNode,*ptptBoundary);
	tmpPtPt = (*ptptNode)->GetPtSelf();
	(*ptptNode)->SetPtSelf((*ptptBoundary)->GetPtSelf());
	(*ptptBoundary)->SetPtSelf(tmpPtPt);
	arrIntCharBoundaries[schar]--;
};

void SR::Hexagon::DecrementCharOfNode(SR::Node** ptptNode, int schar) {
	// See notes on incrementcharofnode
	SR::Node** tmpPtPt;
	// SR::Node* tmpPt;
	SR::Node** ptptBoundary = GetFirstNode()+arrIntCharBoundaries[schar-1];
	swap(*ptptNode,*ptptBoundary);
	tmpPtPt = (*ptptNode)->GetPtSelf();
	(*ptptNode)->SetPtSelf((*ptptBoundary)->GetPtSelf());
	(*ptptBoundary)->SetPtSelf(tmpPtPt);
	arrIntCharBoundaries[schar-1]++;
};

SR::Node** SR::Hexagon::FirstOfChar(int i) {
	if (i==0) return GetFirstNode();
	else return GetFirstNode()+arrIntCharBoundaries[i-1];
};

SR::Node** SR::Hexagon::LastOfChar(int i) {
	// LastOfChar conforms to usual STL approach of being one past the last
	return GetFirstNode()+arrIntCharBoundaries[i];
};

SR::GridHex::GridHex(SR::ParameterSet& p, SR::Hexagon tmphex, SR::DensityField& homes) {

	sizeVecHexagons = p.GetIntValue("intMaxNoHexagons");
	cerr << sizeVecHexagons << "\n";
	sizeVecNodes = p.GetIntValue("intNoNodes");
	vecNodes = new SR::Node[sizeVecNodes];
	vecHexagon = new SR::Hexagon[sizeVecHexagons];
	for (int i=0;i<sizeVecHexagons;++i) vecHexagon[i]=tmphex;
	vecPtNodesHexOrder = new SR::Node*[sizeVecNodes];
	sizeVecHexagons = p.GetIntValue("intMaxNoHexagons");
	cerr << "Entering gridhex constructor...\n";
	maxdx = p.GetValue("dblXGridSize");
	maxdy = p.GetValue("dblYGridSize");
	double minx = p.GetValue("dblXGridMin");
	double miny = p.GetValue("dblYGridMin");
	double dblPoissonAve, dblMaxRelativeDensity, dblAcceptProbability;
	dblPoissonAve = p.GetValue("dblAverageHousehold") - 1;
	SR::Node *tmpPtNode, *tmpPtNodeEndHouse;
	int intHouseholdCounter;
	double tmpx,tmpy;
	tmpPtNode = vecNodes;

	SR::Hexagon *tmpPtHex = vecHexagon;
	SR::Hexagon *ptHexagon;
	SR::Node **tmpPtPtNode;
	SR::Node **ptptNode;

	// cerr << "here" << endl;
	while (tmpPtHex != vecHexagon+sizeVecHexagons) {
		tmpPtHex->ptGrid = this;
		tmpPtHex++;
	}

	// cerr << "here2" << endl;
	dblMaxRelativeDensity = homes.GetMaxVal();
	cerr << "Choosing household locations and sizes using variable household density...\n";
	while (tmpPtNode != vecNodes+sizeVecNodes) { // checking this routine

// NEED to actually pass the RNG variable to the compiler in the make file!
// #ifndef RNG
// #error not defined
// #endif

// #if RNG == gsl
		tmpx = gsl_rng_uniform(glob_rng)*GetMaxDx() + minx;
		tmpy = gsl_rng_uniform(glob_rng)*GetMaxDy() + miny;
// #error one

// #elif RNG == nrcpp

//		tmpx = NR::ran2(p.intSeed)*GetMaxDx() + minx;
//		tmpy = NR::ran2(p.intSeed)*GetMaxDy() + miny;

// #error two

// #else
//		#error Code inbetween NRCPP and GSL random numbers. One must be specified
// #endif

		// cerr << p.intSeed << " " << tmpx << " " << tmpy << " " << dblMaxRelativeDensity << endl;
		dblAcceptProbability = homes.Value(tmpx,tmpy)/dblMaxRelativeDensity;
		// cerr << "here4" << endl;
// 		if (NR::ran2(p.intSeed) < dblAcceptProbability) {
		if (gsl_rng_uniform(glob_rng) < dblAcceptProbability) {
//			intHouseholdCounter = 1 + static_cast<int>(NR::poidev(dblPoissonAve,p.intSeed));
			intHouseholdCounter = 1 + static_cast<int>(gsl_ran_poisson(glob_rng,dblPoissonAve));
			if (vecNodes+sizeVecNodes > tmpPtNode+intHouseholdCounter) tmpPtNodeEndHouse = (tmpPtNode + intHouseholdCounter);
			else {
				tmpPtNodeEndHouse = vecNodes+sizeVecNodes;
				intHouseholdCounter = vecNodes + sizeVecNodes - tmpPtNode;
			}
			while (tmpPtNode != tmpPtNodeEndHouse) {
				tmpPtNode->ptGridHex = this;
				tmpPtNode->SetHouseholdMax(intHouseholdCounter-1);
				tmpPtNode->SetX(tmpx);
				tmpPtNode->SetY(tmpy);
				if (tmpPtNode->GetIndex()%10000==0) cerr << "Nodes assigned: " << tmpPtNode->GetIndex() << "                                  \r";
				tmpPtNode->SetCharacteristic(0);
				tmpPtNode->SetAge(static_cast<int>(gsl_rng_uniform(glob_rng)*80.0));
				tmpPtNode++;
			}
		}
	}
	cerr << "\ndone.\n";

	// Initialise and calculate coordinate extrema
	IntCoord tmpIntCoord;
	intMinXCoord = 9999999;
	intMinYCoord = 9999999;
	dblHexagonWidth = p.GetValue("dblHexagonWidth");
	SR::Node* ptNode = vecNodes;

	tmpIntCoord = RealToHexCoords(GetMaxDx(),GetMaxDy());

	// SR 20140930 these changed line below from +1 to +2. Some approximation in hexagons needs looking at
	// Should this be +!?
	intNoXCoords = tmpIntCoord.x + 2;
	intNoYCoords = tmpIntCoord.y + 2;

	lastHexagon = vecHexagon+intNoXCoords*intNoYCoords;
	cerr << intNoXCoords << "\n";
	cerr << intNoYCoords << "\n";
	if (lastHexagon > vecHexagon+sizeVecHexagons) SR::srerror("Error in calculation of area of hexagon.");

	// Generate all required hexagons (ordering of x,y,node checked by example)
	// Is this section still OK though?
	for (int i=0;i<intNoXCoords;++i) {
		for (int j=0;j<intNoYCoords;++j) {
			(vecHexagon+i*intNoYCoords+j)->SetCoordX(i);
			(vecHexagon+i*intNoYCoords+j)->SetCoordY(j);
		};
	};

	// Alternative method of allocating ptptNodes to hexagons
	// Cycle through each node and increment the right hexagon pointer

	cerr << "Assigning nodes to hexagons and to this gridhex (first node sweep)...";
	ptNode = vecNodes;
	while (ptNode != vecNodes+sizeVecNodes) {
		tmpIntCoord = RealToHexCoords(ptNode->GetX()-minx,ptNode->GetY()-miny);
		(vecHexagon+tmpIntCoord.x*intNoYCoords+tmpIntCoord.y)->IncrementLastNode();
		if (tmpIntCoord.x*intNoYCoords+tmpIntCoord.y >= LastHexagon()-vecHexagon) {
			cerr << "Hexagon allocation out of bounds " << tmpIntCoord.x*intNoYCoords+tmpIntCoord.y << "\n";
			SR::srerror("See previous message");
		}
		ptNode++;
	}
	cerr << "done\n";

	cerr << "Assigning nodes to hexagons (first hexagon sweep)...";
	tmpPtHex = vecHexagon;
	tmpPtPtNode = vecPtNodesHexOrder;

	while (tmpPtHex != LastHexagon()) {
		tmpPtHex->SetFirstNode(tmpPtPtNode);
		tmpPtPtNode += tmpPtHex->ptptLastNode;
		tmpPtHex->SetLastNode(tmpPtHex->GetFirstNode());
		tmpPtHex++;
	}
	cerr << "done\n";

	if (tmpPtPtNode != vecPtNodesHexOrder+sizeVecNodes)
		SR::srerror("Problem allocating nodes to hexagons.");

	cerr << "Assigning nodes to hexagons (second node sweep)...";
	ptNode = vecNodes;
	while (ptNode != vecNodes+sizeVecNodes) {
		tmpIntCoord = RealToHexCoords(ptNode->GetX()-minx,ptNode->GetY()-miny);
		tmpPtHex = vecHexagon+tmpIntCoord.x*intNoYCoords+tmpIntCoord.y;
		tmpPtPtNode = tmpPtHex->GetLastNode();
		*tmpPtPtNode = ptNode;
		tmpPtHex->IncrementLastNode();
		ptNode++;
	}
	cerr << "done\n";

	tmpPtHex = vecHexagon;
	while (tmpPtHex != LastHexagon()-1) {
		if (tmpPtHex->GetLastNode() != (tmpPtHex+1)->GetFirstNode()) {
			cerr << "error";
			exit(1);
		}
		tmpPtHex++;
	};

	// Once memory structures are static (i.e. past here) assign the self pointers and hex pointers
	cerr << "Assigning nodes to hexagons (second hexagon sweep)...";
	ptHexagon = vecHexagon;
	while (ptHexagon != LastHexagon()) {
		ptHexagon->SetLastNode(ptHexagon->GetLastNode()); // Bit strange, but some extras are done in setlastnode
		ptHexagon->AssignCentrePoint(dblHexagonWidth,p);
		ptptNode = ptHexagon->GetFirstNode();
		while (ptptNode!=ptHexagon->GetLastNode()) {
			(*ptptNode)->SetPtSelf(ptptNode);
			(*ptptNode)->ptHexInt = ptHexagon - vecHexagon;
			ptptNode++;
		}
		for (int i=0;i<ptHexagon->usedmn;++i) {
			ptHexagon->SetMaximalNode(vecNodes,i);
		}
		ptHexagon++;
	}
	cerr << "done\n";

	cerr << "Verifying nodes and hexagons...";
	double debugdistance;
	ptNode = vecNodes;
	while (ptNode != vecNodes+sizeVecNodes) {
		debugdistance = ptNode->GetHexagon()->Distance(ptNode);
		if (debugdistance > dblHexagonWidth) {
			// cerr << ptNode - vecNodes << "\t" << debugdistance << "\t" << dblHexagonWidth << "\n";
#ifdef _DEBUG
			// SR::srerror("Problem here with nodes and hexagons.");
#endif
		}
			ptNode++;
	}
	cerr << "done.\n";

	Blocks = 0;
	cerr << "...leaving gridhex constructor.\n";
};

SR::IntCoord SR::GridHex::RealToHexCoords(double xp, double yp) {
	IntCoord rtnpair;
	double xhexoffset=0;
	double yhexoffset=0;
	double x = 1.0 * ( xp - xhexoffset ) / (1.5*dblHexagonWidth);
	double y = 1.0 * ( yp - yhexoffset ) / (dblHexagonWidth*1.7320508075688772935274463415059);
	double z = -0.5 * x - y;
	y = -0.5 * x + y;
	int ix = static_cast<int>(floor(x+0.5));
	int iy = static_cast<int>(floor(y+0.5));
	int iz = static_cast<int>(floor(z+0.5));
	int s = ix+iy+iz;
	if( s )
	{
		double abs_dx = fabs(ix-x);
		double abs_dy = fabs(iy-y);
		double abs_dz = fabs(iz-z);
		if( abs_dx >= abs_dy && abs_dx >= abs_dz )
			ix -= s;
		else if( abs_dy >= abs_dx && abs_dy >= abs_dz )
			iy -= s;
		else
			iz -= s;
	}
	rtnpair.x=ix;
	rtnpair.y=( iy - iz + (1-ix%2) ) / 2;
	return rtnpair;
};

string SR::Hexagon::OutputAllIndexAndLocations() {
	ostringstream oss;
	SR::Node** ptptNode = GetFirstNode();
	while (ptptNode != GetLastNode()) {
		oss << (*ptptNode)->GetIndex() << "\t" << (*ptptNode)->GetX() << "\t" << (*ptptNode)->GetY() << "\n";
		ptptNode++;
	}
	return oss.str();
};

void SR::GridHex::ReserveMemoryForNetwork(SR::PagesForThings<SR::Node>* blks) {
	SR::Node* ptN = vecNodes;
	Blocks=blks;
	int count,msgcount=0;
	cerr << "\n";
	while (ptN!=vecNodes+sizeVecNodes) {
		count =  ptN->GetHouseholdMax();
		count += ptN->GetNoSpatialNeighbour();
		ptN->ptFirstHouseholdMember = Blocks->InsertThing(count);
		ptN++;msgcount++;
		if (msgcount%100000==0) cerr << "Completed " << msgcount << " nodes." << "      \r";
	}
	cerr << "\n";
};

int SR::GridHex::CalculateSizeOfNetwork() {
	SR::Node* ptN = vecNodes;
	int count=0;
	while (ptN!=vecNodes+sizeVecNodes) {
		count +=  ptN->GetHouseholdMax();
		count += ptN->GetNoSpatialNeighbour();
		ptN++;
	}
	return count;
};

void SR::Node::AddToNeighbourList(SR::Node* ptNode1) {
	SR::Node** ptptStart = GetFirstHouseholdMember();
	SR::Node** ptptEnd = ptptStart+GetCurrentPointer();
	while (ptptStart!=ptptEnd) {
		ptptStart++;
	}
	*ptptStart=ptNode1;
	IncrementCurrentPointer();
	return;
};

SR::Node::Node(int i,double x,double y) {
#ifdef SR_BYTEPACKED
	bp1.Set(1,1,0);
	bp1.Set(2,7,0);
	bp1.Set(8,13,0);
	bp1.Set(14,18,0);
	bp2.Set(1,8,0);
	bp2.Set(9,16,0);
	bp2.Set(17,17,0);
	bp2.Set(18,32,0);
#else
	intVaccinationClass=0;
	intCharacteristic=0;
	generation=0;
	no1=0;no2=0;
	intCurrentPointer=0;
	intNoLevelsQuarantine=0;
	blContactsFlag=0;
	fltContactAverage=0;
    intAge=99;
#endif
	dblX=x;
	dblY=y;
};

SR::Node::Node() {
#ifdef SR_BYTEPACKED
	bp1.Set(1,1,0);
	bp1.Set(2,7,0);
	bp1.Set(8,13,0);
	bp1.Set(14,18,0);
	bp2.Set(1,8,0);
	bp2.Set(9,16,0);
	bp2.Set(17,17,0);
	bp2.Set(18,32,0);
#else
	intVaccinationClass=0;
	intCharacteristic=0;
	generation=0;
	no1=0;no2=0;
	intCurrentPointer=0;
	intNoLevelsQuarantine=0;
	blContactsFlag=0;
	fltContactAverage=0;
    intAge=99;
#endif
	dblX=-1;
	dblY=-1;
};

bool SR::GridHex::WriteNodeLocationsToFile(string fn) {
	SR::Node* ptNode = FirstNode();
	SR::Node* ptEnd = LastNode();
	ofstream ofs;
	ofs.open(fn.c_str());
	if (ofs.fail()) return false;
	while (ptNode!=ptEnd) {
		ofs << ptNode->GetIndex() << "\t" << ptNode->GetX() << "\t" << ptNode->GetY() << "\n";
		ptNode++;
	}
	ofs.close();
	return true;
};

bool SR::GridHex::WriteNodeLocationsAndSizesToFile(string fn) {
	SR::Node* ptNode = FirstNode();
	SR::Node* ptEnd = LastNode();
	ofstream ofs;
	ofs.open(fn.c_str());
	int housemembersleft, houseindex=0;
	if (ofs.fail()) return false;

	// Write the header line
	ofs << "index" << ", " << "coord.x" << ", " <<
		"coord.y" << ", " << "age" << ", " << "household.index" << "\n";

	// XXXX this loop is just Wrong
	housemembersleft = ptNode->GetHouseholdMax()+1;
	while (ptNode!=ptEnd) {
		ofs << ptNode->GetIndex() << ", " << ptNode->GetX() << ", " <<
			ptNode->GetY() << ", " << ptNode->GetAge() << ", " << houseindex << "\n";
		ptNode++;
		housemembersleft--;
		if (ptNode != ptEnd && housemembersleft==0) {
			housemembersleft = ptNode->GetHouseholdMax()+1;
			houseindex++;
		}
	}

	ofs.close();
	return true;
};

bool SR::GridHex::WriteArcsToFile(string fn) {
	static int stream_precision=6;
	// SR::Node* tmpNode;
	SR::Node* ptNode = FirstNode();
	SR::Node* ptEnd = LastNode();
	SR::Node **ptptNode,**ptptLast;
	ofstream ofs;
	ofs.open(fn.c_str());
	ofs.precision(stream_precision);
	if (ofs.fail()) return false;

	// Setup the header file
	ofs << "a.index" << ", " << "a.coord.x" << ", " << "a.coord.y" << ", " <<
		"b.index" << ", " << "b.coord.x" << ", " << "b.coord.y" << "\n";

	// Loop through the nodes and non-household links
	while (ptNode != ptEnd) {
		ptptNode = ptNode->GetFirstHouseholdMember()+ptNode->GetHouseholdMax();
		ptptLast = ptptNode+ptNode->GetNoSpatialNeighbour();
		while (ptptNode!=ptptLast) {
			ofs << ptNode->GetIndex() << ", " << ptNode->GetX() << ", " << ptNode->GetY() << ", " <<
				   (*ptptNode)->GetIndex() << ", " << (*ptptNode)->GetY() << ", " << (*ptptNode)->GetY() << "\n";
			ptptNode++;
		}
		ptNode++;
	}

	// Close off the file
	ofs.close();
	return true;
};

bool SR::GridHex::WriteNodeLocationsAndLinksToFile(string fn) {
	// static int stream_precision=6;
	// SR::Node* tmpNode;
	SR::Node* ptNode = FirstNode();
	SR::Node* ptEnd = LastNode();
	ofstream ofs;
	ofs.open(fn.c_str());
	if (ofs.fail()) return false;
	while (ptNode != ptEnd) {
		ofs << ptNode->OutputIndexLocationNetwork() << "\n";
		ptNode++;
	}
	ofs.close();
	return true;
};

SR::Node** SR::Node::GetFirstHouseholdMember() {
	return GetHexagon()->GetPtGrid()->GetPtBlocks()->ReturnThing(ptFirstHouseholdMember);
};

string SR::GridHex::OutputCharacteristicsOnOneLine() {
	ostringstream oss;
	SR::Node* ptNode = FirstNode();
	SR::Node* ptNodeEnd = LastNode();
	while (ptNode != ptNodeEnd) {
		oss << ptNode->GetCharacteristic() << "\t";
		ptNode++;
	}
	oss << "\n";
	return oss.str();
};

string SR::Hexagon::OutputIndexByCharacteristicOnOneLine() {
	ostringstream oss;
	SR::Node** ptptNode = GetFirstNode();
	SR::Node** ptptNodeEnd = GetLastNode();
	oss << "cbs ";
	for (int i=0;i<usedcb;++i) {
		oss << arrIntCharBoundaries[i] << " ";
	}
	oss << "crs ";
	while (ptptNode != ptptNodeEnd) {
		oss << (*ptptNode)->GetIndex() << " " << (*ptptNode)->GetCharacteristic() << " ";
		ptptNode++;
	}
	oss << "\n";
	return oss.str();
};

string SR::GridHex::OutputHexagonCharacteristicsOnManyLines() {
	ostringstream oss;
	SR::Hexagon* ptHexagon = vecHexagon;
	SR::Hexagon* ptHexagonEnd = LastHexagon();
	while (ptHexagon != ptHexagonEnd) {
		oss << ptHexagon->OutputIndexByCharacteristicOnOneLine();
		ptHexagon++;
	}
	return oss.str();
};

void SR::Node::MakeCharacteristicEqualTo(int nc) {
	if (GetCharacteristic() < nc) {
		while (GetCharacteristic() != nc) {
			GetHexagon()->IncrementCharOfNode(GetPtSelf(),GetCharacteristic());
			IncrementBitwiseCharacteristic();
		}
	} else if (GetCharacteristic() > nc) {
		while (GetCharacteristic() != nc) {
			GetHexagon()->DecrementCharOfNode(GetPtSelf(),GetCharacteristic());
			DecrementBitwiseCharacteristic();
		}
	}
};

void SR::GridHex::MakeAllNodesThisCharacteristicAndGeneration(int c, int g) {
/*
	SR::Node* ptNode = vecNodes;
	SR::Node* ptNodeEnd = vecNodes+sizeVecNodes;
	while (ptNode != ptNodeEnd) {
		ptNode->MakeCharacteristicEqualTo(c);
		ptNode->SetGeneration(g);
		ptNode->SetVaccinationClass(0);
		ptNode->ZeroQuarantineLevel();
		ptNode->UnSetContactsFlag();
		ptNode++;
	}
*/
	Hexagon *ptHex = vecHexagon;
	while (ptHex != LastHexagon()) {
		ptHex->ResetMembersCharGen(c,g);
		ptHex++;
	}
};

void SR::GridHex::MakeAllHexagonsNotInRegionalTreatment() {
	static SR::Hexagon* ptHexLocal;
	ptHexLocal = FirstHexagon();
	while (ptHexLocal != LastHexagon()) {
		ptHexLocal->SetRegionalTreatment(0);
		ptHexLocal++;
	}
};

SR::Node* SR::Hexagon::GetMaximalNode(int index) {
#ifdef _DEBUG
	if (index < 0 || index > usedmn) {
		SR::srerror("index out of range in GetMaximalNode.\n");
	}
#endif
	return arrMaximalNodes+index;
};

void SR::Hexagon::SetMaximalNode(SR::Node*  mn, int index) {
#ifdef _DEBUG
	if (index < 0 || index > usedmn) SR::srerror("index out of range in SetMaximalNode.\n");
#endif
	arrMaximalNodes[index] = *mn;
	arrMaximalNodes[index].dblX = GetX();
	arrMaximalNodes[index].dblY = GetY();
	// Node* debug = &arrMaximalNodes[index];
	// int debug2 = 0;
};

double SR::Hexagon::Distance(SR::Node* ptN) {
	return SR::Distance(GetX(),GetY(),ptN->GetX(),ptN->GetY(),GetPtGrid()->GetMaxDx(),GetPtGrid()->GetMaxDy());
	// return sqrt((GetX()-ptN->GetX())*(GetX()-ptN->GetX())+(GetY()-ptN->GetY())*(GetY()-ptN->GetY()));
};

double SR::GridHex::CalculateAverageSpatialNeighbours() {
	SR::Node* ptNode = vecNodes;
	SR::Node* ptNodeEnd = vecNodes+sizeVecNodes;
	double current_total=0;
	while (ptNode != ptNodeEnd) {
		current_total += static_cast<double>(ptNode->GetNoSpatialNeighbour());
		ptNode++;
	}
	current_total /= static_cast<double>(sizeVecNodes);
	return current_total;
};

double SR::GridHex::CalculateAverageHousehold() {
	SR::Node* ptNode = vecNodes;
	SR::Node* ptNodeEnd = vecNodes+sizeVecNodes;
	double current_total=0;
	while (ptNode != ptNodeEnd) {
		current_total += static_cast<double>(ptNode->GetHouseholdMax());
		ptNode++;
	}
	current_total /= static_cast<double>(sizeVecNodes);
	return current_total;
};

string SR::Node::OutputNodeAndHexagonDetails() {
	ostringstream oss;
	oss.precision(6);
	oss << GetIndex() << "\t" << GetX() << "\t" << GetY() << "\t"
		<< GetHexagon() << "\t" << GetHexagon()->GetX() << "\t" << GetHexagon()->GetY() << "\n";
		// << GetHexagon()->Distance(this) << "\n";
	return oss.str();
};

void SR::Hexagon::SetLastNode(SR::Node** ptpt) {
	ptptLastNode=ptpt-GetPtGrid()->GetFirstOfVecPtsHexagon();
	for (int i=0;i<usedcb;++i)
		arrIntCharBoundaries[i]=ptptLastNode-ptptFirstNode;
};

void SR::GridHex::AssignHouseholds() {
	SR::Node *ptCurrentStart,*ptCurrent,*ptCurrentEnd;
	ptCurrentStart=vecNodes;
	// int countmsg=0;
	while (ptCurrentStart != vecNodes+sizeVecNodes) {
		ptCurrentEnd = ptCurrentStart+ptCurrentStart->GetHouseholdMax()+1;
		while (ptCurrentStart != ptCurrentEnd-1) {
			ptCurrent = ptCurrentStart+1;
			while (ptCurrent != ptCurrentEnd) {
				ptCurrentStart->AddToNeighbourList(ptCurrent);
				ptCurrent->AddToNeighbourList(ptCurrentStart);
				ptCurrent++;
			}
			ptCurrentStart++;
		}
		ptCurrentStart=ptCurrentEnd;
	}
};

SR::GridHex::~GridHex() {
	delete [] vecNodes;
	delete [] vecHexagon;
	delete [] vecPtNodesHexOrder;
}




//==================================== READ from/WRITE to binary ===============================




ofstream& SR::operator<<(ofstream& ofs, Node& n) {
	// static char *filePointer;
#ifdef SR_BYTEPACKED
	SR::BinWrite(ofs,n.bp1.GetBasicInt());
	SR::BinWrite(ofs,n.bp2.GetBasicInt());
#else
	SR::BinWrite(ofs,n.intVaccinationClass);
	SR::BinWrite(ofs,n.intCharacteristic);
	SR::BinWrite(ofs,n.generation);
	SR::BinWrite(ofs,n.no1);
	SR::BinWrite(ofs,n.no2);
	SR::BinWrite(ofs,n.intCurrentPointer);
	SR::BinWrite(ofs,n.intNoLevelsQuarantine);
	SR::BinWrite(ofs,n.blContactsFlag);
	SR::BinWrite(ofs,n.fltContactAverage);
    SR::BinWrite(ofs,n.intAge);
	SR::BinWrite(ofs,n.intKernelIndex);
#endif
	SR::BinWrite(ofs,n.dblX);
	SR::BinWrite(ofs,n.dblY);
	SR::BinWrite(ofs,n.ptFirstHouseholdMember.GetPageIndex());
	SR::BinWrite(ofs,n.ptFirstHouseholdMember.GetPtThing());
	SR::BinWrite(ofs,n.ptHexInt);
	SR::BinWrite(ofs,n.ptSelfInt);
	return ofs;
};

ifstream& SR::operator>>(ifstream& ifs, SR::Node& n)
{

	#ifdef SR_BYTEPACKED
		tmp = SR::BinRead<unsigned int>(ifs); n.bp1.SetBasicInt(tmp);
		tmp = SR::BinRead<unsigned int>(ifs); n.bp2.SetBasicInt(tmp);
	#else
		n.intVaccinationClass = SR::BinRead<int>(ifs);
		n.intCharacteristic = SR::BinRead<int>(ifs);
		n.generation = SR::BinRead<int>(ifs);
		n.no1 = SR::BinRead<int>(ifs);
		n.no2 = SR::BinRead<int>(ifs);
		n.intCurrentPointer = SR::BinRead<int>(ifs);
		n.intNoLevelsQuarantine = SR::BinRead<int>(ifs);
		n.blContactsFlag = SR::BinRead<bool>(ifs);
		n.fltContactAverage = SR::BinRead<float>(ifs);
	    n.intAge = SR::BinRead<int>(ifs);
		n.intKernelIndex = SR::BinRead<int>(ifs);
	#endif
		n.dblX = SR::BinRead<float>(ifs);
		n.dblY = SR::BinRead<float>(ifs);
		n.ptFirstHouseholdMember.SetPageIndex(SR::BinRead<int>(ifs));
		n.ptFirstHouseholdMember.SetPtThing(SR::BinRead<int>(ifs));
		n.ptHexInt = SR::BinRead<int>(ifs);
		n.ptSelfInt = SR::BinRead<int>(ifs);

	return ifs;
}


ofstream& SR::operator<<(ofstream& ofs, Hexagon& h) {
	static int nochars;
	static int nomaximals;
	nochars = h.usedcb;
	nomaximals = h.usedmn;
	SR::BinWrite(ofs,nochars);
	SR::BinWrite(ofs,nomaximals);
	SR::BinWrite(ofs,h.blMembersAltered);
	SR::BinWrite(ofs,h.dblCenX);
	SR::BinWrite(ofs,h.dblCenY);
	SR::BinWrite(ofs,h.intCoordX);
	SR::BinWrite(ofs,h.intCoordY);
	SR::BinWrite(ofs,h.ptptFirstNode);
	SR::BinWrite(ofs,h.ptptLastNode);
	SR::BinWrite(ofs,h.intRegionalTreatmentStatus);
	for (int i=0;i<h.usedcb;++i)
	{
		SR::BinWrite(ofs,h.arrIntCharBoundaries[i]);
	}
	for (int i=0;i<h.usedmn;++i)
	{
		ofs << h.arrMaximalNodes[i];
	}

	return ofs;
};

ifstream& SR::operator>>(ifstream& ifs, SR::Hexagon& h)
{
	h.usedcb = SR::BinRead<int>(ifs);
	h.usedmn = SR::BinRead<int>(ifs);
	h.blMembersAltered = SR::BinRead<bool>(ifs);
	h.dblCenX = SR::BinRead<double>(ifs);
	h.dblCenY = SR::BinRead<double>(ifs);
	h.intCoordX = SR::BinRead<int>(ifs);
	h.intCoordY = SR::BinRead<int>(ifs);
	h.ptptFirstNode = SR::BinRead<int>(ifs);
	h.ptptLastNode = SR::BinRead<int>(ifs);
	h.intRegionalTreatmentStatus = SR::BinRead<int>(ifs);

	for (int i=0;i<h.usedcb;++i)
	{
		h.arrIntCharBoundaries[i]=SR::BinRead<int>(ifs);
	}
	for (int i=0;i<h.usedmn;++i)
	{
		//Does this end up duplicating nodes since this cannot see GridHex and is not in vecNodes?
		ifs >> h.arrMaximalNodes[i];
	}
	return ifs;
}


ofstream& SR::operator<<(ofstream& ofs, SR::GridHex& gh) {

	SR::BinWrite(ofs,gh.sizeVecHexagons);
	SR::BinWrite(ofs,gh.sizeVecNodes);
	SR::BinWrite(ofs,gh.maxdx);
	SR::BinWrite(ofs,gh.maxdy);
	SR::BinWrite(ofs,gh.dblHexagonWidth);
	SR::BinWrite(ofs,gh.intMinXCoord);
	SR::BinWrite(ofs,gh.intMinYCoord);
	SR::BinWrite(ofs,gh.intNoXCoords);
	SR::BinWrite(ofs,gh.intNoYCoords);

	for (int i=0;i<gh.sizeVecNodes;++i)
	{
		ofs << gh.vecNodes[i];
	}
	for (int i=0;i<gh.LastHexagon()-gh.FirstHexagon();++i)
	{
		ofs << gh.vecHexagon[i];
	}

	int tmpint;
	for (int i=0;i<gh.sizeVecNodes;++i)
	{
		tmpint = gh.vecPtNodesHexOrder[i]-gh.vecNodes;
		SR::BinWrite(ofs,tmpint);
	}
	return ofs;
};


ifstream& SR::operator>>(ifstream& ifs, SR::GridHex& gh)
{
	gh.sizeVecHexagons = SR::BinRead<int>(ifs);
	gh.sizeVecNodes = SR::BinRead<int>(ifs);
	gh.maxdx = SR::BinRead<double>(ifs);
	gh.maxdy = SR::BinRead<double>(ifs);
	gh.dblHexagonWidth = SR::BinRead<double>(ifs);
	gh.intMinXCoord = SR::BinRead<int>(ifs);
	gh.intMinYCoord = SR::BinRead<int>(ifs);
	gh.intNoXCoords = SR::BinRead<int>(ifs);
	gh.intNoYCoords = SR::BinRead<int>(ifs);

	gh.vecNodes = new SR::Node[gh.sizeVecNodes];
	for (int i=0;i<gh.sizeVecNodes;++i)
	{
		ifs >> gh.vecNodes[i];
		gh.vecNodes[i].SetPtGridHex(&gh);
	}

	gh.vecHexagon = new SR::Hexagon[gh.sizeVecHexagons];
	gh.lastHexagon = gh.vecHexagon+(gh.intNoXCoords*gh.intNoYCoords);
	for (int i=0;i<gh.LastHexagon()-gh.FirstHexagon();++i)
	{
		ifs >> gh.vecHexagon[i];
		gh.vecHexagon[i].SetPtGrid(&gh);
		for (int j=0;j<gh.vecHexagon[i].GetUsedMn();++j)
		{
			//?? Doesn't this get set by hexagon?
			//?? Doesn't seem to return to original value
			//gh.vecHexagon[i].SetMaximalNode(gh.vecNodes,j);
		}
	}

	gh.vecPtNodesHexOrder = new SR::Node*[gh.sizeVecNodes];
	int tmpint;
	for (int i=0;i<gh.sizeVecNodes;++i)
	{
		tmpint = SR::BinRead<int>(ifs);
		gh.vecPtNodesHexOrder[i] = gh.vecNodes+tmpint;
	}

	//This needs reseting?
	gh.Blocks=0;

	return ifs;
}





SR::Node SR::ReadNodeBinaryFromFile(ifstream& ifs) {
	// static char *filePointer;
	// static unsigned int tmp;
	static SR::Node n;
#ifdef SR_BYTEPACKED
	tmp = SR::BinRead<unsigned int>(ifs); n.bp1.SetBasicInt(tmp);
	tmp = SR::BinRead<unsigned int>(ifs); n.bp2.SetBasicInt(tmp);
#else
	n.intVaccinationClass = SR::BinRead<int>(ifs);
	n.intCharacteristic = SR::BinRead<int>(ifs);
	n.generation = SR::BinRead<int>(ifs);
	n.no1 = SR::BinRead<int>(ifs);
	n.no2 = SR::BinRead<int>(ifs);
	n.intCurrentPointer = SR::BinRead<int>(ifs);
	n.intNoLevelsQuarantine = SR::BinRead<int>(ifs);
	n.blContactsFlag = SR::BinRead<bool>(ifs);
	n.fltContactAverage = SR::BinRead<float>(ifs);
        n.intAge = SR::BinRead<int>(ifs);
	n.intKernelIndex = SR::BinRead<int>(ifs);
#endif
	n.dblX = SR::BinRead<float>(ifs);
	n.dblY = SR::BinRead<float>(ifs);
	n.ptFirstHouseholdMember.SetPageIndex(SR::BinRead<int>(ifs));
	n.ptFirstHouseholdMember.SetPtThing(SR::BinRead<int>(ifs));
	n.ptHexInt = SR::BinRead<int>(ifs);
	n.ptSelfInt = SR::BinRead<int>(ifs);
	return n;
};

SR::Hexagon SR::ReadHexagonBinaryFromFile(ifstream& ifs) {
	int tmpNoC, tmpNoM;
	tmpNoC = SR::BinRead<int>(ifs);
	tmpNoM = SR::BinRead<int>(ifs);
	SR::Hexagon h(tmpNoC, tmpNoM);
	h.blMembersAltered = SR::BinRead<bool>(ifs);
	h.dblCenX = SR::BinRead<double>(ifs);
	h.dblCenY = SR::BinRead<double>(ifs);
	h.intCoordX = SR::BinRead<int>(ifs);
	h.intCoordY = SR::BinRead<int>(ifs);
	h.ptptFirstNode = SR::BinRead<int>(ifs);
	h.ptptLastNode = SR::BinRead<int>(ifs);
	h.intRegionalTreatmentStatus = SR::BinRead<int>(ifs);
	for (int i=0;i<h.usedcb;++i) h.arrIntCharBoundaries[i]=SR::BinRead<int>(ifs);
	for (int i=0;i<h.usedmn;++i) h.arrMaximalNodes[i]=SR::ReadNodeBinaryFromFile(ifs);
	return h;
};

SR::GridHex::GridHex(SR::ParameterSet& p, Hexagon tmphex, ifstream& ifs) {
	// sizeVecHexagons = p.GetIntValue("intMaxNoHexagons");
	// sizeVecNodes = p.GetIntValue("intNoNodes");
	sizeVecHexagons = SR::BinRead<int>(ifs);
	sizeVecNodes = SR::BinRead<int>(ifs);
	vecNodes = new SR::Node[sizeVecNodes];
	vecHexagon = new SR::Hexagon[sizeVecHexagons];
	for (int i=0;i<sizeVecHexagons;++i) vecHexagon[i]=tmphex;
	vecPtNodesHexOrder = new SR::Node*[sizeVecNodes];
	SR::Node dbgNode;
	// maxdx = p.GetValue("dblXGridSize");
	// maxdy = p.GetValue("dblYGridSize");
	maxdx = SR::BinRead<double>(ifs);
	maxdy = SR::BinRead<double>(ifs);
	SR::Node* debug;
	dblHexagonWidth = SR::BinRead<double>(ifs);
	intMinXCoord = SR::BinRead<int>(ifs);
	intMinYCoord = SR::BinRead<int>(ifs);
	intNoXCoords = SR::BinRead<int>(ifs);
	intNoYCoords = SR::BinRead<int>(ifs);
	lastHexagon = vecHexagon+intNoXCoords*intNoYCoords;
	SR::Node* dbit;
	int tmpint;
	for (int i=0;i<sizeVecNodes;++i) {
		vecNodes[i] = ReadNodeBinaryFromFile(ifs);
		vecNodes[i].ptGridHex = this;
		dbit = &vecNodes[i];
		dbgNode = vecNodes[i];
	}
	for (int i=0;i<LastHexagon()-FirstHexagon();++i) {
		vecHexagon[i] = ReadHexagonBinaryFromFile(ifs);
		vecHexagon[i].ptGrid = this;
		for (int j=0;j<vecHexagon[i].usedmn;++j) {
			vecHexagon[i].SetMaximalNode(vecNodes,j);
		}
		// cerr << "Hexagons\t" << vecHexagon[i].GetCoordX() << "\t" << vecHexagon[i].GetCoordY() << "\n";
	}
	for (int i=0;i<sizeVecNodes;++i) {
		tmpint = SR::BinRead<int>(ifs);
		vecPtNodesHexOrder[i] = vecNodes+tmpint;
		debug = vecPtNodesHexOrder[i];
		// cerr << "ptsNodes\t" << vecPtNodesHexOrder[i]->GetIndex() << "\n";
	}
	Blocks=0;
};

string SR::Node::OutputNodeToLine() {
	ostringstream oss;
	oss << "\t" << GetIndex() << "\t" << dblX << "\t" << dblY << "\t"
		<< "H\t" << GetHexagon() << "\tS\t" << GetPtSelf() << "\t"
		<< ptFirstHouseholdMember.GetPageIndex() << "\t" << ptFirstHouseholdMember.GetPtThing() << "\n";
	return oss.str();
};

string SR::Hexagon::OutputHexagonToMultipleLines() {
	ostringstream oss;
	oss << "\t" << dblCenX << "\t" << dblCenY << "\t" << intCoordX << "\t" << intCoordY << "\t"
		<< "F\t" << ptptFirstNode << "\tL\t" << ptptLastNode << "\tG\t" << ptGrid << "\t";
	for (int i=0;i<usedcb;++i) oss << arrIntCharBoundaries[i] << "\t";
	oss << "\n";
	oss << "\tMaximal_Nodes:\n";
	for (int i=0;i<usedmn;++i) oss << "\t" << arrMaximalNodes[i].OutputNodeToLine();
	return oss.str();
};

string SR::GridHex::OutputGridHexToMultipleLines() {
	ostringstream oss;
	oss << dblHexagonWidth << "\t" << intMinXCoord << "\t" << intMinYCoord << "\t" << intNoXCoords << "\t" << intNoYCoords << "\n";
	oss << "L\t" << lastHexagon << "\tB\t" << Blocks << "\t" << sizeVecHexagons << "\t" << GetNoHexagons() << "\n";
	oss << "vecNodes:\n";
	for (int i=0;i<sizeVecNodes;++i) oss << "\t" << &vecNodes[i] << vecNodes[i].OutputNodeToLine();
	oss << "vecHexagons:\n";
	for (int i=0;i<sizeVecHexagons;++i) oss << "\t" << &vecHexagon[i] << vecHexagon[i].OutputHexagonToMultipleLines();
	oss << "vecPtsToHexagons:\n";
	for (int i=0;i<sizeVecNodes;++i) oss << "\t" << &vecPtNodesHexOrder[i] << "\t" << vecPtNodesHexOrder[i] << "\t" << vecPtNodesHexOrder[i]->GetIndex() << "\n";
	return oss.str();
};

SR::Node** SR::Hexagon::GetFirstNode() {
	return GetPtGrid()->GetFirstOfVecPtsHexagon()+ptptFirstNode;
};

void SR::Hexagon::SetFirstNode(SR::Node** ptpt) {
	ptptFirstNode=ptpt-GetPtGrid()->GetFirstOfVecPtsHexagon();
};

SR::Node** SR::Hexagon::GetLastNode() {
	return GetPtGrid()->GetFirstOfVecPtsHexagon()+ptptLastNode;
};

int SR::Hexagon::GetNoNodes() {
	return GetLastNode() - GetFirstNode();
};

double SR::Node::Distance(SR::Node* ptN) {
	return SR::Distance(GetX(),GetY(),ptN->GetX(),ptN->GetY(),GetHexagon()->GetPtGrid()->GetMaxDx(),GetHexagon()->GetPtGrid()->GetMaxDy());
};

SR::Hexagon* SR::Node::GetHexagon() {
	return GetPtGridHex()->FirstHexagon()+ptHexInt;
};

void SR::Node::SetPtSelf(SR::Node** pt) {
	ptSelfInt = pt-GetPtGridHex()->GetFirstOfVecPtsHexagon();
};

SR::Node** SR::Node::GetPtSelf() {
	return GetPtGridHex()->GetFirstOfVecPtsHexagon()+ptSelfInt;
};

void SR::GridHex::AssignToPagesOfNodes(SR::PagesForThings<SR::Node>* blks) {
	Blocks = blks;
}

int SR::Node::GetIndex() {
	return this-GetPtGridHex()->FirstNode();
};

void SR::AssignHexagons(SR::Hexagon* ptHex, SR::Hexagon* ptHexEnd, GridHex &g) {
	while (ptHex!=ptHexEnd) {
		// cerr << &g << "\n";
		ptHex->SetPtGrid(&g);
		ptHex++;
	};
};

double SR::GridHex::CalcExpectedSpatial(SR::ParameterSet& p, SR::KERNEL k) {
	// note - currently assumes interactions to be symmetric
	int tmpSeed = -1234;
	cerr << "Calculating expected numbers of infections for spatial kernel...\n";
	// double tmpval = p.GetValue("Relative_Transmit_Spatial");
	p.ChangeValue("Relative_Transmit_Spatial",1);
	double nonodes = static_cast<double>(GetNoNodes());
	int nosamples = static_cast<int>(pow(static_cast<double>(10.0),p.GetIntValue("intPowerTenInteractionsSampled")));
	SR::Node *pt1,*pt2;
	double cumprob=0,rtnval;
	int counter=0;
	for (int i=0;i<nosamples;++i) {
//		pt1 = FirstNode()+static_cast<int>(NR::ran2(tmpSeed)*nonodes);
//		pt2 = FirstNode()+static_cast<int>(NR::ran2(tmpSeed)*nonodes);
		pt1 = FirstNode()+static_cast<int>(gsl_rng_uniform(glob_rng)*nonodes);
		pt2 = FirstNode()+static_cast<int>(gsl_rng_uniform(glob_rng)*nonodes);
		cumprob+=k(p,pt1,pt2,0);
		if (++counter%1000000==0) {
			rtnval = (cumprob/i)*(nonodes-1);
			cerr << "\r" << 100.0*i/static_cast<double>(nosamples) << "% Completed_Sample\t" << " R_0^S* (not scaled) " << rtnval
			<< "                         ";
		}
	}
	cerr << "\ndone.\n";
	rtnval = (cumprob/static_cast<double>(nosamples))*(nonodes-1);
	cerr << "\nNodes:\t" << nonodes << "\tExpected:\t" << rtnval << "\n";
	return rtnval;
};

void SR::Hexagon::ResetMembersCharGen(int c, int g) {
	if (!GetMembersAltered()) return;
	SR::Node** ptptNode = GetFirstNode();
	while (ptptNode != GetLastNode()) {
		(*ptptNode)->MakeCharacteristicEqualTo(c);
		(*ptptNode)->SetGeneration(g);
		(*ptptNode)->SetVaccinationClass(0);
		(*ptptNode)->ZeroQuarantineLevel();
		(*ptptNode)->UnSetContactsFlag();
		ptptNode++;
	}
	UnSetMembersAltered();
}

int SR::Hexagon::GetNoInfectedNotInfectious() {
	static int rtnval;
	rtnval = LastOfChar(3)-FirstOfChar(1)+LastOfChar(8)-FirstOfChar(8)+LastOfChar(7)-FirstOfChar(6);
	return rtnval;
};

int SR::GridHex::GetTotalInfectedNotInfectious() {
	static SR::Hexagon* ptHexLocal;
	static int rtnval;
	ptHexLocal = FirstHexagon();
	rtnval=0;
	while (ptHexLocal != LastHexagon()) {
		rtnval+=ptHexLocal->GetNoInfectedNotInfectious();
		ptHexLocal++;
	}
	return rtnval;
};

SR::NodeMask::NodeMask(SR::GridHex &GH) {

	maxNoNodes = GH.GetNoNodes();
	mask = new bool[maxNoNodes];
	vecNodesSeen = new int[maxNoNodes];
	int seenNoNodes = 0;

};

SR::NodeMask::~NodeMask() {

	delete [] mask;
	delete [] vecNodesSeen;

};


void SR::NodeMask::AgeMask(int lbInc, int ubInc, SR::GridHex &GH) {

	int checkSize, tmpAge;
	SR::Node* ptFirstNode;
	SR::Node* ptTmpNode;
	seenNoNodes = 0;

	// Check for obvious incompatibilities
	checkSize = GH.GetNoNodes();
	if (checkSize != maxNoNodes) SR::srerror("Stopping because mask and gridhex are not compatible");
	ptFirstNode = GH.FirstNode();

	for (int i=0; i<maxNoNodes; ++i) {

		// Check for ages and then add to the bool mask and the lookup mask
		ptTmpNode = ptFirstNode + i;
		tmpAge = ptTmpNode->GetAge();

		if (tmpAge >= lbInc && tmpAge <= ubInc) {
			vecNodesSeen[seenNoNodes] = i;
			seenNoNodes++;
			mask[i] = true;
		} else {
			mask[i] = false;
		}

	}

};

void SR::NodeMask::NullMask(SR::GridHex &GH) {

	// Check for obvious incompatibilities
	int checkSize;
	checkSize = GH.GetNoNodes();
	if (checkSize != maxNoNodes) SR::srerror("Stopping because mask and gridhex are not compatible");

	seenNoNodes = checkSize;

	for (int i=0; i<maxNoNodes; ++i) {
		vecNodesSeen[i] = i;
		mask[i] = true;
	}

};
