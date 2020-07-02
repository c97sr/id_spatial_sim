#include"SR_Workplaces.h"

extern gsl_rng * glob_rng;

void SR::GenerateAllNeighbours(SR::GridHex& gh, SR::ParameterSet& p, PROCESS proc1, KERNEL kern1, KERNEL kern2, UNTIMEDEVENT ue, SR::Workplaces &w) {
	double probability;
	bool eventResult;
	double *ptMSN = p.GetPointer("dblMaxSpatialNeighbours");
	double *ptPCIN = p.GetPointer("dblPropColleguesInNetwork");
	double dblHalfHexOffset = p.GetValue("dblHexagonWidth")/2.0;
	double constIntro = p.GetValue("Constant_Generate_Spatial_Neighbour");
	vector<int> ObjectChar(1,0);
	vector<int> TargetChar(1,0);
	SR::UntimedEvent tmpUe;
	SR::EventMatrix em(1,static_cast<int>(*ptMSN),0,0);
	SR::SpatialKernel spKern(ObjectChar,TargetChar,kern1,proc1,p.GetIntValue("intMaxNeighbourEvents"),p.GetIntValue("intIndexOfMaximalElement"));
	SR::Node *ptNodeCurrent = gh.FirstNode();
	SR::Hexagon* ptHex;
	int sizeCurrentChoices = 1000;
	SR::Node** currentChoices = new SR::Node*[sizeCurrentChoices];
	int currentChoicesUsed=0;
	SR::Node **ptCurrentChoicesFriendA,**ptCurrentChoicesFriendB;
	int *itIntA, *itIntB, *itIntC;
	int ncnt=0; // NMF
	if (constIntro > 0) {
		cerr << "Generating neighbours from spatial kernel...\n";
		while (ptNodeCurrent != gh.LastNode()) {
			ptHex = gh.FirstHexagon();
			em.ClearAllEvents();
			currentChoicesUsed=0;
			while (ptHex != gh.LastHexagon()) {
				spKern.OneNodeOneHex(em,ptHex,ptNodeCurrent,p,dblHalfHexOffset);
				ptHex++;
			}
			em.RemoveDuplicatesAtCurrentTime();
			while (em.CurrentTimeHasEvents()) { //
				tmpUe = em.PopTop();
				currentChoices[currentChoicesUsed]=tmpUe.ptNode2;
				currentChoicesUsed++;
				if (currentChoicesUsed == sizeCurrentChoices) SR::srerror("Current choices problem");
				eventResult = tmpUe.ue(tmpUe.ptNode1,tmpUe.ptNode2,em,p);
			};
			if (currentChoicesUsed>1) {
				ptCurrentChoicesFriendB = currentChoices+1;
				while (ptCurrentChoicesFriendB != currentChoices+currentChoicesUsed) {
					ptCurrentChoicesFriendA = currentChoices;
					while (ptCurrentChoicesFriendA!=ptCurrentChoicesFriendB) {
						probability = kern2(p,*ptCurrentChoicesFriendA,*ptCurrentChoicesFriendB,0);
//						if (NR::ran2(p.intSeed)<probability) eventResult = ue(*ptCurrentChoicesFriendA,*ptCurrentChoicesFriendB,em,p);
						if (gsl_rng_uniform(glob_rng)<probability) eventResult = ue(*ptCurrentChoicesFriendA,*ptCurrentChoicesFriendB,em,p);
						ptCurrentChoicesFriendA++;
					}
					ptCurrentChoicesFriendB++;
				}
			}
			if ((++ncnt)%1000==0) cerr << "Node " << ncnt << "\t" << ptNodeCurrent->GetNoSpatialNeighbour() << "\t"  << "\r"; // NMF
			ptNodeCurrent++;
		}
		cerr << "\ndone.\n";
	}

	// Put p.GetValue("dblPropColleguesInNetwork") proportion of collegues into first level neighbours
	itIntA = w.VeryFirstCollegue();
	ptNodeCurrent = gh.FirstNode();
	cerr << "Generating neighbours from workplaces...\n";
	//DO SOMETHING LESS FRAGILE THAN THIS MAGIC NUMBER
	while (*itIntA > 0 && itIntA != w.OnePastVeryLastCollegue()) {
		itIntB = w.OnePastLastCollegue(*itIntA);
		itIntC = itIntA+1;
		while (itIntC != itIntB) {
//			if (NR::ran2(p.intSeed)<*ptPCIN) {
//			XXXX for some reason, if I swap in the line below, I start to get seg faults. No idea why.
			if ((static_cast<double>(gsl_rng_uniform(glob_rng))) < (*ptPCIN)) {
				eventResult = ue(ptNodeCurrent+*itIntA,ptNodeCurrent+*itIntC,em,p);
			}
			itIntC++;
		}
		itIntA++;
		if ((itIntA-w.VeryFirstCollegue()) % 100000==0) cerr << "Node " << itIntA-w.VeryFirstCollegue() << "           \r";
	}

	delete [] currentChoices;

	cerr << "\ndone.\n";

};

void SR::Workplaces::MakeWorkplacesConsistentWithNodes(SR::GridHex& g, SR::ParameterSet& p) {
	int currentNode,tmpint;
	NodeMask *nm = g.GetMaskPtr();
	SR::workplace* ptFirstWorkplace = vecWorkplaces;
	SR::workplace* ptOnePastLastWorkplace = vecWorkplaces+sizeVecWorkplaces;
	SR::workplace* ptWorkplace;
	SR::Node* ptFirstHexNode = g.FirstNode();
	SR::Node* ptOnePastLastHexNode = g.LastNode();
	SR::Node* ptHexNode;
	// Allocate space by setting counter ints
	cerr << "\nMaking workplaces consistent with nodes. Starting sweep through workplaces to set pointers...";
	ptWorkplace = ptFirstWorkplace;
	tmpint = 0;
	while (ptWorkplace!=ptOnePastLastWorkplace) {
		ptWorkplace->firstnode = tmpint;
		tmpint += ptWorkplace->totalnodes;
		ptWorkplace->currentlyusednodes=0;
		ptWorkplace++;
	}
	cerr << "  done.\n";

	// Repeat node pass assigning the node this time
	cerr << "  Starting second sweep through nodes assigning workplaces...";

	for(int i=0; i<nm->GetNoSeenNodes(); i++) {
		ptHexNode = g.FirstNode() + nm->RevealNode(i);
		currentNode = ptHexNode->GetIndex();
		tmpint = vecWorkplaceIntsForEachNode[currentNode];
		ptWorkplace = ptFirstWorkplace+tmpint;
		vecNodeIntsInWorkplaceOrder[ptWorkplace->firstnode+ptWorkplace->currentlyusednodes]=currentNode;
		ptWorkplace->currentlyusednodes++;
		ptHexNode++;
	}

}

SR::Workplaces::~Workplaces() {
	delete [] vecNodeIntsInWorkplaceOrder;
	delete [] vecWorkplaceIntsForEachNode;
	delete [] vecWorkplaces;
	delete [] vecDistFrequencies;
	delete [] vecDistFrequenciesOld;
	delete [] ptFirstOrderedWorkplace;
	delete [] ptIntFirstNoWorkplacesInGrid;
	delete [] ptFirstPtWorkPlace;
	delete [] ptArrNumberNearbyWorkplaces;
	delete [] ptArrPtPtNearbyWorkplaces;
	delete [] ptArrNumberEachNearbyWorkplaces;
}

SR::Workplaces::Workplaces(SR::GridHex& g, SR::ParameterSet& p, SR::DensityField& wpd) {

	cerr << "Entering workplaces contructor...\n";

	// Old inheritance bits
	sizeVecNodeIntsInWorkplaceOrder = p.GetIntValue("intNoNodes");
	sizeVecWorkplaceIntsForEachNode = p.GetIntValue("intNoNodes");
	sizeVecWorkplaces = p.GetIntValue("intNumberOfWorkplaces");
	sizeVecDistFrequencies = outputBinNumber;

	vecNodeIntsInWorkplaceOrder = new int[sizeVecNodeIntsInWorkplaceOrder];
	vecWorkplaceIntsForEachNode = new int[sizeVecNodeIntsInWorkplaceOrder];
	for (int i=0;i<sizeVecNodeIntsInWorkplaceOrder;++i) vecNodeIntsInWorkplaceOrder[i]=-1;
	for (int i=0;i<sizeVecWorkplaceIntsForEachNode;++i) vecWorkplaceIntsForEachNode[i]=-1;

	vecWorkplaces = new SR::workplace[sizeVecWorkplaces];
	vecDistFrequencies = new int[sizeVecDistFrequencies];
	vecDistFrequenciesOld = new int[sizeVecDistFrequencies];
	for (int i=0;i<sizeVecDistFrequencies;++i) {vecDistFrequencies[i]=0;vecDistFrequenciesOld[i]=1;}
	degreesfreedom = static_cast<int>(sqrt(g.GetMaxDx()*g.GetMaxDx()+g.GetMaxDy()*g.GetMaxDy()));


	outputBinSize=1.0;
	startscale=0.001;
	dblGridSize = p.GetValue("dblMCMCGridSize");
	intNoNodes = g.GetNoNodes();
	intNoWorkplaces = p.GetIntValue("intNumberOfWorkplaces");
	intNoXWorkplaceGrids = static_cast<int>(p.GetValue("dblXGridSize")/dblGridSize+1);
	intNoYWorkplaceGrids = static_cast<int>(p.GetValue("dblYGridSize")/dblGridSize+1);
	intNoWorkPlaceGrids = intNoXWorkplaceGrids*intNoYWorkplaceGrids;
	ptFirstOrderedWorkplace = new SR::workplace*[intNoWorkplaces];
	ptIntFirstNoWorkplacesInGrid = new int[intNoWorkPlaceGrids]; for (int i=0;i<intNoWorkPlaceGrids;++i) ptIntFirstNoWorkplacesInGrid[i]=0;
	ptFirstPtWorkPlace = new SR::workplace**[intNoWorkPlaceGrids];
	int* tmpWorkPlaceGrid = new int[intNoWorkPlaceGrids];
	ptArrNumberNearbyWorkplaces = new int[intNoWorkPlaceGrids]; for (int i=0;i<intNoWorkPlaceGrids;++i) ptArrNumberNearbyWorkplaces[i]=0;
	ptArrPtPtNearbyWorkplaces = new SR::workplace**[intNoWorkPlaceGrids*9]; for (int i=0;i<9*intNoWorkPlaceGrids;++i) ptArrPtPtNearbyWorkplaces[i]=0;
	ptArrNumberEachNearbyWorkplaces = new int[intNoWorkPlaceGrids*9]; for (int i=0;i<9*intNoWorkPlaceGrids;++i) ptArrNumberEachNearbyWorkplaces[i]=0;
	int tmpXint, tmpYint;

	// Initialise required variables
	double xmin = p.GetValue("dblXGridMin");
	double xmax = p.GetValue("dblXGridMin")+p.GetValue("dblXGridSize");
	double ymin = p.GetValue("dblYGridMin");
	double ymax = p.GetValue("dblYGridMin")+p.GetValue("dblYGridSize");
	double dblDenseMax = wpd.GetMaxVal();
	double dblDenseCurrent, dblRelProb, tmpx,tmpy, tmpdist;
	int intSeedLog,tmpint,intSquare2,intSquare1;
	SR::workplace* ptFirstWorkplace = vecWorkplaces;
	SR::workplace* ptOnePastLastWorkplace = vecWorkplaces+sizeVecWorkplaces;
	SR::workplace* ptWorkplace;
	SR::Node* ptFirstHexNode = g.FirstNode();
	SR::Node* ptOnePastLastHexNode = g.LastNode();
	SR::Node* ptHexNode;
	ptWorkplace = ptFirstWorkplace;
	int indexcounter=0;

	// Initialise all workplaces with an x and y coords. Set counters.  Count numbers of workplaces in each grid space.
	cerr << "Assigning spatial locations to " << p.GetIntValue("intNumberOfWorkplaces") << " workplaces...\n";
	while (ptWorkplace!=ptOnePastLastWorkplace) {
//		tmpx = NR::ran2(p.intSeed)*(xmax-xmin)+xmin;
//		tmpy = NR::ran2(p.intSeed)*(ymax-ymin)+ymin;
		tmpx = gsl_rng_uniform(glob_rng)*(xmax-xmin)+xmin;
		tmpy = gsl_rng_uniform(glob_rng)*(ymax-ymin)+ymin;
		dblDenseCurrent = wpd.Value(tmpx,tmpy);
		dblRelProb = dblDenseCurrent / dblDenseMax;
//		if (NR::ran2(p.intSeed) <= dblRelProb) {
		if (gsl_rng_uniform(glob_rng) <= dblRelProb) {
			ptWorkplace->x = tmpx;
			ptWorkplace->y = tmpy;
			ptWorkplace->totalnodes = 0;
			ptWorkplace->currentlyusednodes=0;
			ptWorkplace->index = indexcounter;
			indexcounter++;
			ptWorkplace++;
			tmpXint = static_cast<int>((tmpx-xmin)/dblGridSize);
			tmpYint = static_cast<int>((tmpy-ymin)/dblGridSize);
			ptIntFirstNoWorkplacesInGrid[tmpXint*intNoYWorkplaceGrids+tmpYint]++;
			if (indexcounter%100000==0) cerr << "Workplaces assigned: " << indexcounter << "                                  \r";
		}
	}

	tmpint = 0;
	for (int i=0;i<intNoWorkPlaceGrids;++i) {
		tmpWorkPlaceGrid[i] = tmpint;
		ptFirstPtWorkPlace[i] = ptFirstOrderedWorkplace+tmpint;
		tmpint += ptIntFirstNoWorkplacesInGrid[i];
	}
	if (tmpint != intNoWorkplaces) SR::srerror("Problem in allocation of workplace pointers.");

	indexcounter = 0;
	cerr << "\n";
	ptWorkplace =  vecWorkplaces;
	while (ptWorkplace!=ptOnePastLastWorkplace) {
		tmpx = ptWorkplace->x;
		tmpy = ptWorkplace->y;
		tmpXint = static_cast<int>((tmpx-xmin)/dblGridSize);
		tmpYint = static_cast<int>((tmpy-ymin)/dblGridSize);
		tmpint = tmpXint*intNoYWorkplaceGrids + tmpYint;
		ptFirstOrderedWorkplace[tmpWorkPlaceGrid[tmpint]]=ptWorkplace;
		tmpWorkPlaceGrid[tmpint]++;
		indexcounter++;
		if (indexcounter%100000==0) cerr << "Workplaces assigned to grids " << indexcounter << "                                  \r";
		ptWorkplace++;
	}

	cerr << "\n";

	for (int i=0;i<intNoWorkPlaceGrids;++i) tmpWorkPlaceGrid[i]=0;

	for (int i=0;i<intNoWorkPlaceGrids;++i) {
		intSquare1 = i;
		intSquare2 = intSquare1 - intNoYWorkplaceGrids;
		while (intSquare2 != intSquare1 + 2*intNoYWorkplaceGrids) {
			if (intSquare2 >= 0 && intSquare2 < intNoWorkPlaceGrids) {
				ptArrNumberNearbyWorkplaces[i] += ptIntFirstNoWorkplacesInGrid[intSquare2];
				ptArrNumberEachNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]] = ptIntFirstNoWorkplacesInGrid[intSquare2];
				ptArrPtPtNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]]=ptFirstPtWorkPlace[intSquare2];
				tmpWorkPlaceGrid[i]++;
				// update list of pointers to pointers
				// cerr << intNoYWorkplaceGrids << "\t" << intSquare1  << "\t" << intSquare2 << "\n";
			}
			intSquare2 += intNoYWorkplaceGrids; // this either has to be inside or out, not both
		}

		if (i%intNoYWorkplaceGrids!=0) {
			intSquare2 = intSquare1 - 1 - intNoYWorkplaceGrids;
			while (intSquare2 != intSquare1 -1 + 2*intNoYWorkplaceGrids) {
				if (intSquare2 >= 0 && intSquare2 < intNoWorkPlaceGrids) {
					ptArrNumberNearbyWorkplaces[i] += ptIntFirstNoWorkplacesInGrid[intSquare2];
					ptArrNumberEachNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]] = ptIntFirstNoWorkplacesInGrid[intSquare2];
					ptArrPtPtNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]]=ptFirstPtWorkPlace[intSquare2];
					tmpWorkPlaceGrid[i]++;
					// update list of pointers to pointers
					// cerr << intNoYWorkplaceGrids << "\t" << intSquare1  << "\t" << intSquare2 << "\n";
				}
				intSquare2 += intNoYWorkplaceGrids; // this either has to be inside or out, not both
			}
		}

		if (i%intNoYWorkplaceGrids!=intNoYWorkplaceGrids-1) {
			intSquare2 = intSquare1 + 1 - intNoYWorkplaceGrids;
			while (intSquare2 != intSquare1 + 1 + 2*intNoYWorkplaceGrids) {
				if (intSquare2 >= 0 && intSquare2 < intNoWorkPlaceGrids) {
					ptArrNumberNearbyWorkplaces[i] += ptIntFirstNoWorkplacesInGrid[intSquare2];
					ptArrNumberEachNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]] = ptIntFirstNoWorkplacesInGrid[intSquare2];
					ptArrPtPtNearbyWorkplaces[i*9+tmpWorkPlaceGrid[i]]=ptFirstPtWorkPlace[intSquare2];
					tmpWorkPlaceGrid[i]++;
					// update list of pointers to pointers
					// cerr << intNoYWorkplaceGrids << "\t" << intSquare1  << "\t" << intSquare2 << "\n";
				}
				intSquare2 += intNoYWorkplaceGrids; // this either has to be inside or out, not both
			}
		}

	}

	NodeMask *nm = g.GetMaskPtr();
	for(int i=0; i<nm->GetNoSeenNodes(); i++) {
		//cerr << nm->RevealNode(i);
		ptHexNode = g.FirstNode() + nm->RevealNode(i);
		tmpint = static_cast<int>(gsl_rng_uniform(glob_rng)*static_cast<double>(sizeVecWorkplaces));
		vecWorkplaceIntsForEachNode[ptHexNode->GetIndex()]=tmpint;


		ptWorkplace = ptFirstWorkplace+tmpint;
		tmpdist = SR::Distance(ptHexNode->GetX(),ptHexNode->GetY(),ptWorkplace->x,ptWorkplace->y,g.GetMaxDx(),g.GetMaxDy());
		tmpint = GetDistanceBin(tmpdist);
		vecDistFrequencies[tmpint]++;
		(ptWorkplace->totalnodes)++;
	}
	cerr << "  ...done.\n";

	MakeWorkplacesConsistentWithNodes(g, p);

	delete [] tmpWorkPlaceGrid;

	cerr << "done\n... leaving workplaces contructor.\n";


	// DEBUG
	/*tmpint = 0;
	int sum = 0;
	cerr <<"\n\n\n";
	cerr << "STATE OF vecNodeIntsInWorkplaceOrder";
	for(int i = 0; i < sizeVecNodeIntsInWorkplaceOrder; i++)
	{
		if (i >= sum) {
			cerr << "\nWP" << tmpint << " Nodes:\n";
			sum += (ptFirstWorkplace+tmpint)->totalnodes;
			tmpint+=1;
		}
		cerr << vecNodeIntsInWorkplaceOrder[i] <<" ";
	}
	cerr <<"\n\n\n";
	cerr << "STATE OF vecWorkplaceIntsForEachNode\n";
	for(int i = 0; i < sizeVecNodeIntsInWorkplaceOrder; i++)
	{
		cerr << "Node:" << i << " WP:" << vecWorkplaceIntsForEachNode[i] << "\n";
	}

*/

};

SR::workplace* SR::Workplaces::GetWorkplaceOfNode(int ni) {
	return vecWorkplaces+vecWorkplaceIntsForEachNode[ni];
	static int tmpint;
	if (ni < 0 || ni >= sizeVecWorkplaceIntsForEachNode) SR::srerror("Problem in SR::Workplaces::GetWorkplaceOfNode(int ni)");
	tmpint = vecWorkplaceIntsForEachNode[ni];
	if (tmpint < 0 || tmpint >= sizeVecWorkplaces) {
		SR::srerror("Problem in SR::Workplaces::GetWorkplaceOfNode(int ni) 2nd prob");
	}
	SR::workplace* rtnit = vecWorkplaces+tmpint;
	return rtnit;
};

void SR::Workplaces::WriteWorkplacesToFile(string f) {
	ofstream ofs;
	SR::workplace *it = vecWorkplaces;
	ofs.open(f.c_str());
	if (ofs.fail()) SR::srerror("Problem in SR::Workplaces::WriteWorkplacesToFile(string f).");
	while (it!= vecWorkplaces+sizeVecWorkplaces) {
		ofs << it->index << "\t" << it->x << "\t" << it->y << "\t" << it->totalnodes << "\n";
		it++;
	}
	ofs.close();
};

void SR::Workplaces::WriteCommutesToFile(SR::GridHex& g, string f) {
	ofstream ofs;
	SR::workplace* itwk;
	SR::Node* itnd = g.FirstNode();
	double dist;
	ofs.open(f.c_str());
	if (ofs.fail()) SR::srerror("Problem in SR::Workplaces::WriteWorkplacesToFile(string f).");
	while (itnd != g.LastNode()) {
		itwk = GetWorkplaceOfNode(itnd->GetIndex());
		dist = SR::Distance(itnd->GetX(),itnd->GetY(),itwk->x,itwk->y, g.GetMaxDx(), g.GetMaxDy());
		ofs << itnd->GetX() << "\t" << itnd->GetY() << "\t" << itwk->x << "\t" << itwk->y << "\t" << dist << "\n";
		itnd++;
	}
	ofs.close();
};

void SR::Workplaces::MCMCUpdate(
		GridHex& gh,
		ParameterSet& p,
		double pdistance(double,int,double*),
		string funcfile,
		int minSamplesMillions,
		int maxSamplesMillions,
		int& cn) {

	NodeMask *nm = gh.GetMaskPtr();
	bool blMoveAccepted;
	double dblAcceptanceProbability,dblProposalProbability;
	int intSourceBinIndex,intTargetBinIndex;
	int *ptIntSourceBin,*ptIntTargetBin;
	double dblCurrentCommute, dblProposedCommute, dblP_i, dblP_j;
	bool notYetStable = true;
	int indexSelectedNode, indexMaskedNode, indexProposedWorkplace, intCount=0;
	SR::Node *ptSelectedNode;
	SR::workplace *ptCurrentWorkplace, *ptProposedWorkplace;
	double aveDeltaRatio = 0;
	double aveDeltaRatioSquared = 0;
	double epsilon;
	double dblWeightLocal = p.GetValue("dblWeightMCMCUpdateLocal");
	double dblLocalThreshold = dblWeightLocal / (dblWeightLocal + 1);
	double xoffset = p.GetValue("dblXGridMin");
	double yoffset = p.GetValue("dblYGridMin");
	pair<double,double> chiTest;
	SR::CachedDblLookup *ptLookup;

	int intNoAllWorkplaces = sizeVecWorkplaces;
	int intNoNearbyWorkplaces;

	int noNodesSeen = nm->GetNoSeenNodes();

	// These might be the wrong indexes... 0 and 1, no?
	double vecp0[2];
	vecp0[0] = p.GetValue("Commute_Power_One_GZ");
	vecp0[1] = p.GetValue("Commute_Change_Point_GZ");

	double vecp1[2];
	vecp1[0] = p.GetValue("Commute_Power_One_HK");
	vecp1[1] = p.GetValue("Commute_Change_Point_HK");
	SR::CachedDblLookup cached_pdistance_0(pdistance,startscale,incsPerOrder,noOrders,2,vecp0);
	SR::CachedDblLookup cached_pdistance_1(pdistance,startscale,incsPerOrder,noOrders,2,vecp1);

	// SR::OpenNullFile(funcfile,cached_pdistance.PrintDistributionOfCommutes());

	double movecount=0;
	int tmpIndexWorkplace;
	int intTmpNodeGridIndex;
	int intKingsSquareOffset;
	pair<int,int> tmpNode,tmpCurrentWorkplace,tmpTargetWorkplace;
	bool blOnlyTargetClose=false,blOnlyCurrentClose=false,blTargetClose,blCurrentClose;

	cerr << "Starting MCMC updates to network...\n";

	while ((notYetStable && intCount < maxSamplesMillions) || (intCount < minSamplesMillions)) {
		for (int i=0;i<1000000;++i) {
			blMoveAccepted=false;

			// select one node and one (different) workplace
//			indexSelectedNode = static_cast<int>(1.0*gh.GetNoNodes()*NR::ran2(p.intSeed));
                        // XXXX Up to here 9th May 2020. Think here is one place I need to alter
                        // and then the place they get assigned initially.
			indexMaskedNode = static_cast<int>(1.0*noNodesSeen*gsl_rng_uniform(glob_rng));
			indexSelectedNode = nm->RevealNode(indexMaskedNode);
			ptSelectedNode = gh.FirstNode()+indexSelectedNode;


			ptCurrentWorkplace = vecWorkplaces+vecWorkplaceIntsForEachNode[indexSelectedNode];
			ptProposedWorkplace = ptCurrentWorkplace;

			// Calculate grid x,y
			tmpNode = CalcGridXY(ptSelectedNode->GetX()-xoffset,ptSelectedNode->GetY()-yoffset);
			tmpCurrentWorkplace = CalcGridXY(ptCurrentWorkplace->x-xoffset,ptCurrentWorkplace->y-yoffset);
			intTmpNodeGridIndex = GetGridIndex(tmpNode);
			intNoNearbyWorkplaces = ptArrNumberNearbyWorkplaces[intTmpNodeGridIndex];

			// Why is line -1 "> 0" and line -3 "!=1"
//			if (intNoNearbyWorkplaces > 0 && NR::ran2(p.intSeed) < dblLocalThreshold) {
			if (intNoNearbyWorkplaces > 0 && gsl_rng_uniform(glob_rng) < dblLocalThreshold) {
				// Choose workplace number
				while (ptProposedWorkplace == ptCurrentWorkplace && intNoNearbyWorkplaces != 1) {
					intKingsSquareOffset=0;
//					tmpIndexWorkplace = static_cast<int>(static_cast<double>(intNoNearbyWorkplaces)*NR::ran2(p.intSeed));
					tmpIndexWorkplace = static_cast<int>(static_cast<double>(intNoNearbyWorkplaces)*gsl_rng_uniform(glob_rng));
					while (tmpIndexWorkplace >= ptArrNumberEachNearbyWorkplaces[intTmpNodeGridIndex*9+intKingsSquareOffset]) {
						tmpIndexWorkplace -= ptArrNumberEachNearbyWorkplaces[intTmpNodeGridIndex*9+intKingsSquareOffset];
						intKingsSquareOffset++;
					}
					ptProposedWorkplace = *(ptArrPtPtNearbyWorkplaces[intTmpNodeGridIndex*9+intKingsSquareOffset]+tmpIndexWorkplace);
				}
				indexProposedWorkplace = ptProposedWorkplace->index;
				if (indexProposedWorkplace < 0 || indexProposedWorkplace >= sizeVecWorkplaces) SR::srerror("problem");
			} else {
				while (ptProposedWorkplace==ptCurrentWorkplace) {

					// intNoAllWorkplaces seems badly wrong
//					indexProposedWorkplace = static_cast<int>(1.0*intNoAllWorkplaces*NR::ran2(p.intSeed));
					indexProposedWorkplace = static_cast<int>(1.0*intNoAllWorkplaces*gsl_rng_uniform(glob_rng));
					ptProposedWorkplace = vecWorkplaces+indexProposedWorkplace;
				}
			}

			tmpTargetWorkplace = CalcGridXY(ptProposedWorkplace->x-xoffset,ptProposedWorkplace->y-yoffset);

			blTargetClose = GridIsClose(tmpTargetWorkplace,tmpNode);
			blCurrentClose = GridIsClose(tmpCurrentWorkplace,tmpNode);
			blOnlyTargetClose = blTargetClose && ! blCurrentClose;
			blOnlyCurrentClose = blCurrentClose && ! blTargetClose;

			if (blOnlyTargetClose) {
				dblProposalProbability = 1 / (1+dblLocalThreshold*static_cast<double>(intNoAllWorkplaces) /
					((1-dblLocalThreshold)*static_cast<double>(intNoNearbyWorkplaces)));
			} else if (blOnlyCurrentClose) {
				dblProposalProbability = (1+dblLocalThreshold*static_cast<double>(intNoAllWorkplaces) /
					((1-dblLocalThreshold)*static_cast<double>(intNoNearbyWorkplaces)));
			} else {
				dblProposalProbability = 1;
			}

			// calculate the values required for the acceptance probability
			if (ptSelectedNode->GetKernelIndex()==1) ptLookup = &cached_pdistance_0;
			else ptLookup = &cached_pdistance_1;
			dblCurrentCommute = SR::Distance(ptSelectedNode->GetX(),ptSelectedNode->GetY(),ptCurrentWorkplace->x,ptCurrentWorkplace->y, gh.GetMaxDx(), gh.GetMaxDy());
			dblProposedCommute = SR::Distance(ptSelectedNode->GetX(),ptSelectedNode->GetY(),ptProposedWorkplace->x,ptProposedWorkplace->y, gh.GetMaxDx(), gh.GetMaxDy());
			intSourceBinIndex=GetDistanceBin(dblCurrentCommute);
			intTargetBinIndex=GetDistanceBin(dblProposedCommute);
			ptIntSourceBin=&vecDistFrequencies[intSourceBinIndex];
			ptIntTargetBin=&vecDistFrequencies[intTargetBinIndex];
			dblP_j=ptLookup->GetCheckedValue(dblProposedCommute);
			dblP_i=ptLookup->GetCheckedValue(dblCurrentCommute);
			dblAcceptanceProbability = dblProposalProbability*dblP_j/dblP_i;

			if (dblAcceptanceProbability < 1e-100 || dblAcceptanceProbability > 1e100) {
				SR::srerror("Fallen over because of MCMC acceptance probability");
			}

			if (dblAcceptanceProbability >= 1.0) {blMoveAccepted = true;movecount++;}
//			else if (NR::ran2(p.intSeed) < dblAcceptanceProbability) {blMoveAccepted = true;movecount++;}
			else if (gsl_rng_uniform(glob_rng) < dblAcceptanceProbability) {blMoveAccepted = true;movecount++;}

			// update if move accepted
			if (blMoveAccepted) {
				MoveNode(indexSelectedNode,indexProposedWorkplace);
				aveDeltaRatio+=log(dblP_j/dblP_i);
				aveDeltaRatioSquared+=log(dblP_j/dblP_i)*log(dblP_j/dblP_i);
				(*ptIntSourceBin)--;
				(*ptIntTargetBin)++;
			}

		}


		intCount++;
		cn++;
		aveDeltaRatio /= 1000000; aveDeltaRatioSquared /= 1000000; movecount /= 1000000;
		epsilon = aveDeltaRatio/(sqrt(aveDeltaRatioSquared-aveDeltaRatio*aveDeltaRatio)/1000);
		// chiTest = CalcChiSquared(1.96,degreesfreedom);
		// ArchiveCommuteDist();
		cerr << "\rMU(pdates)\t" << cn << "\tChi_cur\t" << chiTest.first << "\tChi_crt\t" << chiTest.second << "          ";
		// if (chiTest.first < chiTest.second) notYetStable = false;
		aveDeltaRatio = 0; movecount = 0; aveDeltaRatioSquared=0;
	}
	MakeWorkplacesConsistentWithNodes(gh,p);
	cerr << "\n...done.\n";
};

void SR::Workplaces::MoveNode(int nd, int wkf) {
	// static SR::workpalce* currentWorkplaceOfNode, nextWorkplaceOfNode;
	static SR::workplace *currentWorkplaceOfNode, *nextWorkplaceOfNode;
	static int *ptInt;
	ptInt = vecWorkplaceIntsForEachNode+nd;
	currentWorkplaceOfNode = vecWorkplaces+*ptInt;
	nextWorkplaceOfNode = vecWorkplaces+wkf;
	currentWorkplaceOfNode->totalnodes--;
	nextWorkplaceOfNode->totalnodes++;
	// vecWorkplaceIntsForEachNode[nd]=wkf;
	*ptInt=wkf;
}

void SR::Workplaces::MoveNodeDash(int nd, int wkf) {
	static int *currentPosOfNode, *nextPosOfNode;
	static SR::workplace *currentWorkplaceOfNode;
	if (nd >= sizeVecWorkplaceIntsForEachNode || wkf >= sizeVecWorkplaces) SR::srerror("Problem in SR::Workplaces::MoveNode(int nd, wkf)");
	currentWorkplaceOfNode=GetWorkplaceOfNode(nd);
	vecWorkplaceIntsForEachNode[nd]=wkf;
	if (currentWorkplaceOfNode->index < wkf) {
		currentPosOfNode = vecNodeIntsInWorkplaceOrder+currentWorkplaceOfNode->firstnode;
		while (*currentPosOfNode!=nd) currentPosOfNode++;
		nextPosOfNode = vecNodeIntsInWorkplaceOrder+
			currentWorkplaceOfNode->firstnode+
			currentWorkplaceOfNode->totalnodes-1;
		swap(*currentPosOfNode,*nextPosOfNode);
		currentPosOfNode=nextPosOfNode;
		currentWorkplaceOfNode->totalnodes--;
		currentWorkplaceOfNode++;
		currentWorkplaceOfNode->totalnodes++;
		currentWorkplaceOfNode->firstnode--;
		while (currentWorkplaceOfNode->index != wkf) {
			nextPosOfNode = vecNodeIntsInWorkplaceOrder+
				currentWorkplaceOfNode->firstnode+
				currentWorkplaceOfNode->totalnodes-1;
			swap(*currentPosOfNode,*nextPosOfNode);
			currentPosOfNode=nextPosOfNode;
			currentWorkplaceOfNode->totalnodes--;
			currentWorkplaceOfNode++;
			currentWorkplaceOfNode->totalnodes++;
			currentWorkplaceOfNode->firstnode--;
		}
	} else if (currentWorkplaceOfNode->index > wkf) {
		currentPosOfNode = vecNodeIntsInWorkplaceOrder+currentWorkplaceOfNode->firstnode;
		while (*currentPosOfNode!=nd) currentPosOfNode++;
		nextPosOfNode = vecNodeIntsInWorkplaceOrder+currentWorkplaceOfNode->firstnode;
		swap(*currentPosOfNode,*nextPosOfNode);
		currentPosOfNode=nextPosOfNode;
		currentWorkplaceOfNode->totalnodes--;
		currentWorkplaceOfNode->firstnode++;
		currentWorkplaceOfNode--;
		currentWorkplaceOfNode->totalnodes++;
		while (currentWorkplaceOfNode->index != wkf) {
			nextPosOfNode = vecNodeIntsInWorkplaceOrder+currentWorkplaceOfNode->firstnode;
			swap(*currentPosOfNode,*nextPosOfNode);
			currentPosOfNode=nextPosOfNode;
			currentWorkplaceOfNode->totalnodes--;
			currentWorkplaceOfNode->firstnode++;
			currentWorkplaceOfNode--;
			currentWorkplaceOfNode->totalnodes++;
		}
	} else SR::srerror("Trying to move node to itself.");
};

int* SR::Workplaces::FirstCollegue(int n) {
	SR::workplace* wkpl = GetWorkplaceOfNode(n);
	return vecNodeIntsInWorkplaceOrder+wkpl->firstnode;
};

int* SR::Workplaces::OnePastLastCollegue(int n) {
	SR::workplace* wkpl = GetWorkplaceOfNode(n);
	return vecNodeIntsInWorkplaceOrder+wkpl->firstnode+wkpl->totalnodes;
};

int SR::Workplaces::GetDistanceBin(double d) {
	/* Old log scale distance bins
	static double scaleval;
	static double logscalestart = log10(startscale);
	static double dblIncs = static_cast<double>(incsPerOrder);
	static int max = noOrders * incsPerOrder;
	int rtnval;
	scaleval = log10(d) - logscalestart;
	rtnval = static_cast<int>(scaleval*dblIncs)-1;
	if (rtnval < 0) return 0;
	if (rtnval >= max) return max;
	return rtnval;
	*/

	int rtnval;
	rtnval = static_cast<int>(d/outputBinSize);
	if (rtnval < 0) SR::srerror("Negative distance not possible.");
	if (rtnval >= outputBinNumber) rtnval = outputBinNumber-1;
	return rtnval;

};

string SR::Workplaces::PrintDistributionOfCommutes() {
	ostringstream oss;
	oss << "dist.lb, freq\n";
	for (int i=0;i<outputBinNumber;++i) {
		// oss << pow(10,log10(startscale)+i/static_cast<double>(incsPerOrder)) << "\t" << vecDistFrequencies[i] << "\n";
		oss << static_cast<double>(i)*outputBinSize << ", " << vecDistFrequencies[i] << "\n";
	}
	return oss.str();
};

string SR::Workplaces::PrintDistributionOfSizes() {
	ostringstream oss;
	int* vecFrequencyOfSizes;
	workplace* ptWp;
	workplace* ptLastWp = vecWorkplaces + sizeVecWorkplaces;
	int maxsize = 0;
	ptWp = vecWorkplaces;
	while (ptWp != ptLastWp) {
		if (ptWp->totalnodes > maxsize) maxsize = ptWp->totalnodes;
		ptWp++;
	}
	vecFrequencyOfSizes = new int[maxsize+1];
	for (int i=0;i<=maxsize;++i) vecFrequencyOfSizes[i]=0;
	ptWp = vecWorkplaces;
	while (ptWp != ptLastWp) {
		vecFrequencyOfSizes[ptWp->totalnodes]++;
		ptWp++;
	}
	for (int i=0;i<=maxsize;++i) {
		oss << vecFrequencyOfSizes[i] << "\n";
	}
	delete [] vecFrequencyOfSizes;
	return oss.str();
};

pair<int,int> SR::Workplaces::CalcGridXY(double x, double y) {
	pair<int,int> rtnpair;
	rtnpair.first = static_cast<int>(x/dblGridSize);
	rtnpair.second = static_cast<int>(y/dblGridSize);
	return rtnpair;
};

bool SR::Workplaces::GridIsClose(pair<int,int>& a,pair<int,int>& b) {
	if ((abs(a.first-b.first) <= 1) && (abs(a.second-b.second) <= 1)) return true;
	else return false;
};

int SR::Workplaces::GetGridIndex(pair<int,int>& p) {
	return p.first*intNoYWorkplaceGrids + p.second;
};

void SR::Workplaces::GenerateDistributionOfPossibleCommutes(GridHex& gh, ParameterSet& p, double grid_dx, double hist_dx, string funcfile) {

	// decide the dimensions for the grid
	double xmin = p.GetValue("dblXGridMin");
	double ymin = p.GetValue("dblYGridMin");
	int noXCoords = static_cast<int>(p.GetValue("dblXGridSize")/grid_dx+1);
	int noYCoords = static_cast<int>(p.GetValue("dblYGridSize")/grid_dx+1);
	int noHistBins = static_cast<int>(sqrt(static_cast<double>(noXCoords*noXCoords+noYCoords*noYCoords))*grid_dx/hist_dx+1);
	ostringstream oss;

	// Declare the arrays for numbers of people and numbers of workplaces
	vector<double> vecWorkPlaceDist(noHistBins,0);
	vector<int> vecGridWorkplaces(noXCoords*noYCoords,0);
	vector<int> vecGridPopulation(noXCoords*noYCoords,0);

	// Set up the pointers
	SR::workplace* ptWp = vecWorkplaces;
	SR::Node* ptNd = gh.FirstNode();
	int tmpXCoord,tmpYCoord;
	int tmpHistBin;
	double tmpDistance;

	// Set up the grids
	while (ptNd != gh.LastNode()) {
		tmpXCoord = static_cast<int>((ptNd->GetX()-xmin)/grid_dx);
		tmpYCoord = static_cast<int>((ptNd->GetY()-ymin)/grid_dx);
#ifdef _DEBUG
		if (tmpXCoord*noYCoords+tmpYCoord>vecGridPopulation.size()) SR::srerror("Problem with the gridding of workplaces");
#endif
		vecGridWorkplaces[tmpXCoord*noYCoords+tmpYCoord]++;
		ptNd++;
	}
	while (ptWp != vecWorkplaces+sizeVecWorkplaces) {
		tmpXCoord = static_cast<int>((ptWp->x-xmin)/grid_dx);
		tmpYCoord = static_cast<int>((ptWp->y-ymin)/grid_dx);
#ifdef _DEBUG
		if (tmpXCoord*noYCoords+tmpYCoord>vecGridPopulation.size()) SR::srerror("Problem with the gridding of workplaces");
#endif
		vecGridPopulation[tmpXCoord*noYCoords+tmpYCoord]++;
		ptWp++;
	}

	// Calculate the histograms
	cerr << "\n";
	for (int wpXCoord = 0; wpXCoord < noXCoords; wpXCoord++) {
		for (int wpYCoord = 0; wpYCoord < noYCoords; wpYCoord++) {
			for (int ndXCoord = 0; ndXCoord < noXCoords; ndXCoord++) {
				for (int ndYCoord = 0; ndYCoord < noYCoords; ndYCoord++) {
					tmpDistance = sqrt(static_cast<double>(
						(wpXCoord-ndXCoord)*(wpXCoord-ndXCoord) + (wpYCoord-ndYCoord)*(wpYCoord-ndYCoord)
						)) * grid_dx;
					tmpHistBin = static_cast<int>(tmpDistance/hist_dx);
#ifdef _DEBUG
					if (tmpHistBin > vecWorkPlaceDist.size()) SR::srerror("Problem with the histogram part of the this bit.");
#endif
					vecWorkPlaceDist[tmpHistBin]+= static_cast<double>(vecGridPopulation[ndXCoord*noYCoords+ndYCoord]) *
						static_cast<double>(vecGridWorkplaces[wpXCoord*noYCoords+wpYCoord]);
				}
			}
		}
		cerr << wpXCoord+1 << " of " << noXCoords << " completed        \r";
	}

	oss << "dist.lb, freq\n";

	for (int i=0;i<static_cast<int>(vecWorkPlaceDist.size());i++) {
		oss << static_cast<double>(i)*hist_dx << ", " << vecWorkPlaceDist[i] << "\n";
	}

	SR::OpenNullFile(funcfile,oss.str());

};

pair<double, double> SR::Workplaces::CalcChiSquared(double zed_val, int df) {
	int oldResidual=0, newResidual=0;
	int binCount=0;
	double chisquared=0;
	int *ptIntNew,*ptIntOld;
	ptIntNew=vecDistFrequencies;
	ptIntOld=vecDistFrequenciesOld;
	pair<double,double> rtnp(0,0);
	for (int i=0;i<sizeVecDistFrequencies;++i) {
		if (*ptIntOld >= 5) {
			chisquared += static_cast<double>((*ptIntNew-*ptIntOld)*(*ptIntNew-*ptIntOld)) / static_cast<double>(*ptIntOld);
			binCount++;
		} else {
			oldResidual+=*ptIntOld;
			newResidual+=*ptIntNew;
		}
		ptIntNew++;ptIntOld++;
	}
	if (oldResidual >= 5) {
		chisquared += static_cast<double>((newResidual-oldResidual)*(newResidual-oldResidual)) / static_cast<double>(oldResidual);
		binCount++;
	}
	rtnp.first = chisquared;
	rtnp.second = 0.5 * (zed_val+sqrt(2*static_cast<double>(df)-1))*(zed_val+sqrt(2*static_cast<double>(df)-1));
	return rtnp;
};

void SR::Workplaces::ArchiveCommuteDist() {
	for (int i=0;i<sizeVecDistFrequencies;++i) vecDistFrequenciesOld[i]=vecDistFrequencies[i];
}
