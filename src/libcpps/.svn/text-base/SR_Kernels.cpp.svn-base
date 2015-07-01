#include"SR_Kernels.h"

void SR::SpatialKernel::EagerOneNodeOneHex(SR::EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p) {
	double dblCurrentProbability;
	vector<int>::iterator ptTargetChar;
	SR::Node **ptptNodeCurrent, **ptptNodeEnd;
	ptTargetChar=vecIntTargetCharacteristics.begin();
	while (ptTargetChar!=vecIntTargetCharacteristics.end()) {
		ptptNodeCurrent = ptHex->FirstOfChar(*ptTargetChar);
		ptptNodeEnd = ptHex->LastOfChar(*ptTargetChar);
		while (ptptNodeCurrent != ptptNodeEnd) {
			dblCurrentProbability = kernKernel(p,ptN,*ptptNodeCurrent,0);
			if (NR::ran2(p.intSeed)<dblCurrentProbability) {
				procProcess(p,ptN,*ptptNodeCurrent,em);
			}
			ptptNodeCurrent++;
		}
		ptTargetChar++;
	}
};

SR::Kernel::Kernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int mi) :
vecIntSourceCharacteristics(sc.size()), vecIntTargetCharacteristics(sc.size()) {
	vecIntSourceCharacteristics=sc;
	vecIntTargetCharacteristics=tc;
	kernKernel=k;
	procProcess=p;
	intMaximalIndex = mi;
	dblExpNumber = -1.0;
	dblCurrentConst = -2.0;
};

SR::SpatialKernel::SpatialKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int ss, int mi) :
SR::Kernel(sc,tc,k,p,mi){
	sizeofstack = ss;
	vecStackOfInts = new int[sizeofstack];
	intLevelOfStack=0;
};

SR::SpatialKernel::~SpatialKernel(){
	delete [] vecStackOfInts;
};

SR::HouseholdKernel::HouseholdKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int mi) :
SR::Kernel(sc,tc,k,p,mi) {
};

SR::NeighbourKernel::NeighbourKernel(vector<int>& sc, vector<int>& tc, SR::KERNEL k, SR::PROCESS p, int mi) :
SR::Kernel(sc,tc,k,p,mi) {
};

void SR::Kernel::GenerateAllEvents(SR::ParameterSet& p, SR::EventMatrix& em, GridHex& gh) {
	// Changes to this method be propagated to "With Cross" branch.
	SR::Node **ptptNode,**ptptNodeEnd;
	SR::Hexagon *ptHexagonSource;
	SR::Hexagon *ptHexagonEnd = gh.LastHexagon();
	ptHexagonSource = gh.FirstHexagon();
	while (ptHexagonSource != ptHexagonEnd) {
		for (unsigned int i=0;i < vecIntSourceCharacteristics.size();++i) {
			ptptNode = ptHexagonSource->FirstOfChar(vecIntSourceCharacteristics[i]);
			ptptNodeEnd = ptHexagonSource->LastOfChar(vecIntSourceCharacteristics[i]);
			while (ptptNode != ptptNodeEnd) {
				GenerateEventsOneSourceNode(p,em,gh,ptptNode);
				ptptNode++;
			}
		}
		ptHexagonSource++;
	}
};

void SR::Kernel::GenerateAllEventsWithCross(SR::ParameterSet& p, SR::EventMatrix& em, GridHex& gh, SR::ONE_D_DISTRIBUTION c, double* ptC) {
	SR::Node **ptptNode,**ptptNodeEnd;
	SR::Hexagon *ptHexagonSource;
	SR::Hexagon *ptHexagonEnd = gh.LastHexagon();
	ptHexagonSource = gh.FirstHexagon();
	while (ptHexagonSource != ptHexagonEnd) {
		for (unsigned int i=0;i<vecIntSourceCharacteristics.size();++i) {
			ptptNode = ptHexagonSource->FirstOfChar(vecIntSourceCharacteristics[i]);
			ptptNodeEnd = ptHexagonSource->LastOfChar(vecIntSourceCharacteristics[i]);
			while (ptptNode != ptptNodeEnd) {
				GenerateEventsOneSourceNodeWithCross(p,em,gh,ptptNode,c,ptC);
				ptptNode++;
			}
		}
		ptHexagonSource++;
	}
};

void SR::Kernel::GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix& em, GridHex& gh, SR::Node** ptptNode) {
	SR::srerror("Generate events one target node shouldn't have been called.\n");
};

void SR::Kernel::GenerateEventsOneSourceNodeWithCross(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode,SR::ONE_D_DISTRIBUTION condDist, double *ptC) {
	static int intMaxInfsPerDt = p.GetIntValue("intMaxInfsPerDt");
	static int intMaxDelayInAuxEventMatrix = p.GetIntValue("intMaxDelayInAuxEventMatrix");
	static int Maximum_Number_CT_Vaccinations_Each_Day = p.GetIntValue("Maximum_Number_CT_Vaccinations_Each_Day"); // Cannot change between runs at the moment
	static double *dblDailySpatial;
	static double h_rrf = p.GetValue("Hazard_Rash_Relative_Fever");
	static SR::UNTIMEDEVENT nullUe=0;
	static SR::EventMatrix tmpEventMatrix(intMaxDelayInAuxEventMatrix,intMaxInfsPerDt,nullUe,Maximum_Number_CT_Vaccinations_Each_Day);
	static double dblOvercallFactor=2;
	static int intRequiredEvents;
	static double tmpDbl;

	dblDailySpatial = p.GetPointer("Daily_Expected_Spatial");

	if (*dblDailySpatial > 0) {

		if (!(*ptptNode)->GetContactsFlag()) {
			(*ptptNode)->SetContactsFlag();
			tmpDbl = condDist(*dblDailySpatial,p);
			if (tmpDbl < 0) SR::srerror("Problem with conditional distribution.");
			if (tmpDbl > 163) tmpDbl = 163; // greatest int less than 2^(32-18)/100
			(*ptptNode)->SetNoContacts(tmpDbl);
		}

		if ((*ptptNode)->GetCharacteristic()==2) tmpDbl = (*ptptNode)->GetNoContacts();
		else if ((*ptptNode)->GetCharacteristic()==3) tmpDbl = (*ptptNode)->GetNoContacts()*h_rrf;
		else {
			SR::srerror("Node must be in one of these characteristics.");
		}

		intRequiredEvents = static_cast<int>(NR::poidev(tmpDbl,p.intSeed));

		if (intRequiredEvents > 0) {
			*ptC *= intRequiredEvents*dblOvercallFactor;
			while (tmpEventMatrix.TotalEventsPending() < intRequiredEvents) {
				GenerateEventsOneSourceNode(p,tmpEventMatrix,gh,ptptNode);
			}
			*ptC /= intRequiredEvents*dblOvercallFactor;
			while (tmpEventMatrix.TotalEventsPending() != intRequiredEvents) {
				tmpEventMatrix.RemoveEventAtRandom(p);
			}
			tmpEventMatrix.ExportEvents(em);
		}
	}
};

void SR::SpatialKernel::GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix& em, GridHex& gh, SR::Node** ptptNode) {
	SR::Hexagon *ptHexagonTarget = gh.FirstHexagon();
	SR::Hexagon *ptHexagonEnd = gh.LastHexagon();
	double maxRate = p.GetValue("dblMaxWithinHexSampleRate");
	while (ptHexagonTarget != ptHexagonEnd) {
		OneNodeOneHex(em,ptHexagonTarget,*ptptNode,p,maxRate);
		ptHexagonTarget++;
	}
};

void SR::SpatialKernel::LazyOneNodeOneHex(SR::EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double sp) {
	static int intTotalPoss,intNumberChosen,intNodeToTry;
	intLevelOfStack=0;
	intTotalPoss=0;
	SR::Node **ptptStartNode,**ptptCurrentNode;
	ptptStartNode = ptHex->FirstOfChar(0);

	for (unsigned int i=0;i<vecIntTargetCharacteristics.size();++i) {
		intTotalPoss += ptHex->LastOfChar(vecIntTargetCharacteristics[i]) - ptHex->FirstOfChar(vecIntTargetCharacteristics[i]);
	}
	intNumberChosen = static_cast<int>(ignbin(sp,intTotalPoss,p.intSeed));
	if (intNumberChosen > sizeofstack-1) {
		SR::srerror("Too many events chosen in LazyOneNodeOneHex.\n");
	}
	while (intLevelOfStack!=intNumberChosen) {
		intNodeToTry = static_cast<int>(NR::ran2(p.intSeed)*intTotalPoss);
		AddToStack(intNodeToTry);
	}
	if(intNumberChosen>1) sort(vecStackOfInts,vecStackOfInts+intLevelOfStack);
	ptInt = vecStackOfInts;
	ptIntEnd = ptInt+intLevelOfStack;
	if(intNumberChosen>0)
		for (unsigned int i=0;i<vecIntTargetCharacteristics.size();++i) {
			while (ptHex->FirstOfChar(vecIntTargetCharacteristics[i])+*ptInt < ptHex->LastOfChar(vecIntTargetCharacteristics[i]) &&
				ptInt != ptIntEnd) {
					*ptInt += static_cast<int>(ptHex->FirstOfChar(vecIntTargetCharacteristics[i])-ptptStartNode);
					ptInt++;
				}
		}
		ptInt = vecStackOfInts;

		static double debugdistance;
		static double debugdistance2;

		while (ptInt != ptIntEnd) {
			ptptCurrentNode = ptptStartNode + *ptInt;
			probability = kernKernel(p,ptN,*ptptCurrentNode,0);
			if (probability > sp) {
				debugdistance = (*ptptCurrentNode)->Distance(ptN);
				debugdistance2 = (*ptptCurrentNode)->GetHexagon()->Distance(ptN);
				cerr << "\n";
				cerr << (*ptptCurrentNode)->GetX() << "\n";
				cerr << (*ptptCurrentNode)->GetY() << "\n";
				cerr << (*ptptCurrentNode)->GetHexagon()->GetX() << "\n";
				cerr << (*ptptCurrentNode)->GetHexagon()->GetY() << "\n";
				cerr << "This was a problem with sample probability.\n";
				exit(1);
			}
			probability /= sp;
			// double debugdistance = (*ptptCurrentNode)->Distance(ptN);
			// double debugdistance2 = (*ptptCurrentNode)->GetHexagon()->Distance(ptN);
			if (NR::ran2(p.intSeed) < probability) procProcess(p,ptN,*ptptCurrentNode,em);
			ptInt++;
		}
};

void SR::SpatialKernel::AddToStack(int i) {
	ptStack = vecStackOfInts;
	ptStackTop = ptStack + intLevelOfStack;
	while (ptStack != ptStackTop && *ptStack != i) ptStack++;
	if (ptStack==ptStackTop) {*ptStack = i;intLevelOfStack++;}
};

void SR::SpatialKernel::OneNodeOneHex(SR::EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double cp) {
	static double maxprob,initialDistance;
	initialDistance =ptHex->Distance(ptN);
	// if(initialDistance<*ptMaxNodeHexDist) {
	// cout << ptN->GetX() << "\t" << ptN->GetY() << "\t" << ptHex->GetMaximalNode(0)->GetX()
	// 	<< "\t" << ptHex->GetMaximalNode(0)->GetY() << "\t" << initialDistance << "\t" << ptHex->GetPtGrid()->GetHexagonWidth() << "\n";
	if (initialDistance < ptHex->GetPtGrid()->GetHexagonWidth()) {
		EagerOneNodeOneHex(em,ptHex,ptN,p);
	} else {
		maxprob = kernKernel(p,ptN,ptHex->GetMaximalNode(intMaximalIndex),ptHex->GetPtGrid()->GetHexagonWidth());
		if (maxprob > cp) {
			EagerOneNodeOneHex(em,ptHex,ptN,p);
		} else {
			LazyOneNodeOneHex(em,ptHex,ptN,p,maxprob);
		}
	}
	// }
};

void SR::Kernel::OneNodeOneHex(SR::EventMatrix &em, SR::Hexagon *ptHex, SR::Node* ptN, SR::ParameterSet &p, double cp) {
	SR::srerror("Shouldn't have been called.\n");
};

void SR::HouseholdKernel::GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode) {
	Node** ptptHousemate = (*ptptNode)->GetFirstHouseholdMember();
	Node** ptptEndOfHousemates = ptptHousemate + (*ptptNode)->GetHouseholdMax();
	static double probability;
	while (ptptHousemate != ptptEndOfHousemates) {
		probability = kernKernel(p,*ptptNode,*ptptHousemate,0);
		if (NR::ran2(p.intSeed) < probability) procProcess(p,*ptptNode,*ptptHousemate,em);
		ptptHousemate++;
	}
};

void SR::NeighbourKernel::GenerateEventsOneSourceNode(SR::ParameterSet& p, SR::EventMatrix&em, GridHex& gh, SR::Node** ptptNode) {
	Node** ptptHousemate = (*ptptNode)->GetFirstHouseholdMember()+(*ptptNode)->GetHouseholdMax();
	Node** ptptEndOfHousemates = ptptHousemate + (*ptptNode)->GetNoSpatialNeighbour();
	static double probability;
	while (ptptHousemate != ptptEndOfHousemates) {
		probability = kernKernel(p,*ptptNode,*ptptHousemate,0);
		if (NR::ran2(p.intSeed) < probability) procProcess(p,*ptptNode,*ptptHousemate,em);
		ptptHousemate++;
	}
};

SR::HexKernel::HexKernel(vector<int>& sc, PROCESS proc, int mh, double gmr, int gmn, double lmr) :
vecIntSourceCharacteristics(sc.size(),-1) {
	maxHex = mh;
	vecPtHex = new SR::Hexagon*[maxHex];
	vecIntSourceCharacteristics = sc;
	ptOnePastLast = vecPtHex;
	procProcess = proc;
	intPopSize=0;
	globalMaxNoTreatments=gmn;
	globalMaxTreatmentRate=gmr;
	localMaximumRate=lmr;
};

SR::HexKernel::~HexKernel() {
	delete [] vecPtHex;
};

void SR::HexKernel::GenerateEventsForOneHex(SR::Hexagon *ptHex, SR::EventMatrix &em, SR::ParameterSet &p) {
	// Generate vector of pts to pointers of possible nodes
	static int intDosesAvailable;
	static int hexagonPop;
	static bool blPeopleLeft;
	static SR::Node **currentPtPtNode,**KStateEndPtPtNode;
	static int sourceVecCounter;
	static double locallyConstrainedRate,globallyConstrainedRate;

	hexagonPop = ptHex->GetNoNodes();
	locallyConstrainedRate = static_cast<double>(hexagonPop)*localMaximumRate;
	globallyConstrainedRate = static_cast<double>(hexagonPop)/static_cast<double>(intPopSize)*globalMaxTreatmentRate;

	if (locallyConstrainedRate < globallyConstrainedRate) intDosesAvailable = static_cast<int>(locallyConstrainedRate);
	else intDosesAvailable = static_cast<int>(globallyConstrainedRate);

	blPeopleLeft=true;
	sourceVecCounter=0;

	currentPtPtNode=ptHex->FirstOfChar(vecIntSourceCharacteristics[sourceVecCounter]);
	KStateEndPtPtNode=ptHex->LastOfChar(vecIntSourceCharacteristics[sourceVecCounter]);

	while (intDosesAvailable > 0 && blPeopleLeft) {
		if (currentPtPtNode != KStateEndPtPtNode) {
			// Vaccinate or treat
			procProcess(p,*currentPtPtNode,*currentPtPtNode,em);
			intDosesAvailable--;
			currentPtPtNode++;
		} else {
			// Switch to next K State
			sourceVecCounter++;
			if (sourceVecCounter==static_cast<int>(vecIntSourceCharacteristics.size())) {
				blPeopleLeft=false;
			} else {
				currentPtPtNode=ptHex->FirstOfChar(vecIntSourceCharacteristics[sourceVecCounter]);
				KStateEndPtPtNode=ptHex->LastOfChar(vecIntSourceCharacteristics[sourceVecCounter]);
			}
		}
	}
};

void SR::HexKernel::GenerateEvents(EventMatrix &em, ParameterSet& p) {
	SR::Hexagon **ptptHex = vecPtHex;
	double dblSizeCurrentHex;
	while (ptptHex != ptOnePastLast) {
		dblSizeCurrentHex = static_cast<double>((*ptptHex)->GetLastNode()-(*ptptHex)->GetFirstNode());
		GenerateEventsForOneHex(*ptptHex,em,p);
		ptptHex++;
	}
};

void SR::HexKernel::AddHexagons(SR::GridHex &gh, double x, double y, double range) {
	SR::Hexagon *ptHex = gh.FirstHexagon();
	SR::Hexagon *ptLastHex = gh.LastHexagon();
	SR::Hexagon **ptptHex;
	double dblDist;
	while (ptHex != ptLastHex) {
		dblDist = SR::Distance(ptHex->GetX(),ptHex->GetY(),x,y,1e100,1e100);
		if (dblDist <= range && ptHex->GetRegionalTreatment()==0) {
			ptptHex = vecPtHex;
			while (*ptptHex != ptHex && ptptHex != ptOnePastLast) ptptHex++;
			if (ptptHex==ptOnePastLast) {
				intPopSize+=ptHex->GetNoNodes();
				*ptOnePastLast = ptHex;
				ptOnePastLast++;
				if (ptOnePastLast == vecPtHex+maxHex) SR::srerror("Not enough hexagons reserved in HexKernel");
			}
		}
		ptHex++;
	}
};

void SR::HexKernel::RemoveHexagon(SR::Hexagon *ptHex) {
	static SR::Hexagon **ptptHex;
	static bool notyetdone;
	ptptHex = vecPtHex;
	notyetdone=true;
	while (notyetdone && (ptptHex != ptOnePastLast)) {
		if (*ptptHex==ptHex) {
			ptOnePastLast--;
			*ptptHex=*ptOnePastLast;
			notyetdone=false;
			ptHex->SetRegionalTreatment(100);
		}
		ptptHex++;
	}
	if (notyetdone) SR::srerror("Tried to remove a hexagon that wasn't actually in the treatment zone.");
};

void SR::HexKernel::AddHexagons(SR::GridHex &gh, SR::Hexagon *ptHexSource, double range) {
	// Regional treatment in two levels units 0-not in, 1-active and 2-done
	// Tens - 0 not triggered 10+ has triggered
	// This routine zeros out if the source has already been used as a trigger
	static SR::Hexagon *ptHex;
	static SR::Hexagon *ptLastHex;
//	static SR::Hexagon **ptptHex;
//	static int intdebug,intdebug2;
	if (ptHexSource->GetRegionalTreatment()<10) {
		ptHex = gh.FirstHexagon();
		ptLastHex = gh.LastHexagon();
		static double dblDist;
		while (ptHex != ptLastHex) {
			dblDist = SR::Distance(ptHex->GetX(),ptHex->GetY(),ptHexSource->GetX(),ptHexSource->GetY(),1e100,1e100);
			if (dblDist <= range && (ptHex->GetRegionalTreatment()==0 || ptHex->GetRegionalTreatment()==10) && ptHex->GetNoNodes() > 0) {
				if (ptOnePastLast == vecPtHex+maxHex) SR::srerror("Not enough hexagons reserved in HexKernel");
				intPopSize+=ptHex->GetNoNodes();
				*ptOnePastLast = ptHex;
				ptOnePastLast++;
				ptHex->SetRegionalTreatment(ptHex->GetRegionalTreatment()+1);
			}
			ptHex++;
		}
		ptHexSource->SetRegionalTreatment(ptHexSource->GetRegionalTreatment()+10);
	}
};

void SR::HexKernel::UpdateRegion(GridHex &gh, double range, double max_symp, double max_sus) {
	// Hexagons added if they are within
	static SR::Hexagon *ptHex;
	static double prev;
//	static int dblDbg;
	ptHex = gh.FirstHexagon();
	while (ptHex != gh.LastHexagon()) {
		if (ptHex->GetRegionalTreatment() < 10) {
			prev = static_cast<double>(ptHex->LastOfChar(3)-ptHex->FirstOfChar(3)+ptHex->LastOfChar(7)-ptHex->FirstOfChar(7))/
				 static_cast<double>(ptHex->GetNoNodes()+1);
			if (prev > max_symp) {
				AddHexagons(gh,ptHex,range);
			}
		} else if (ptHex->GetRegionalTreatment() < 100) {
			prev = static_cast<double>(ptHex->LastOfChar(0)-ptHex->FirstOfChar(0))/
				static_cast<double>(ptHex->GetNoNodes()+1);
			if (prev < max_sus) {
				RemoveHexagon(ptHex);
			}
		}
		ptHex++;
	}
};

void SR::HexKernel::ResetRegion(GridHex &gh) {
	static SR::Hexagon *ptHex;
	ptHex = gh.FirstHexagon();
	while (ptHex != gh.LastHexagon()) {
		ptHex->SetRegionalTreatment(0);
		ptHex++;
	}
	ptOnePastLast = vecPtHex;
	intPopSize=0;
};

