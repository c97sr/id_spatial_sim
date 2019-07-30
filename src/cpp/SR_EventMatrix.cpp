#include"SR_EventMatrix.h"

extern gsl_rng * glob_rng;

// TODO here . Possible to have delay of -1
void SR::EventMatrix::AddEvent(SR::UntimedEvent ue, int delay) {
	static int initialDelay;
	if (delay >= intMaxDelay) {
		delay = intMaxDelay - 1;
		cerr << "Warning: delay too long\n";
	}
	delay += intCurrentMatrixTime;
	delay  = delay % intMaxDelay;
	initialDelay=delay;
	int flag_full=false; //NMF
	if (ue.ue==constrainedEv) {
		while ((vecNumberOfConstrainedEvents[delay] >= intLevelConstraint)&&(!flag_full)) {
			delay++;
			if (delay==intMaxDelay) delay=0;
			if (delay==initialDelay) flag_full=true;
		}
	}
	if(!flag_full) {
		if (vecSizeOfStacks[delay]==intMaxEvents) {
			cout << vecSizeOfStacks[delay] << "\n";
			SR::srerror("Too many events generated in AddEvents.\n");
		}
		matStacksLinear[delay*intMaxEvents+vecSizeOfStacks[delay]]=ue;
		vecSizeOfStacks[delay]++;intTotalEvents++;
		if (ue.ue==constrainedEv) vecNumberOfConstrainedEvents[delay]++;
	}
};

SR::EventMatrix::EventMatrix(int md, int me, UNTIMEDEVENT cev, int levcon) :
matStacksLinear(md*me),vecSizeOfStacks(md,0),vecNumberOfConstrainedEvents(md,0) {
	intTotalEvents=0;
	intMaxDelay=md;
	intMaxEvents=me;
	intCurrentMatrixTime=0;
	constrainedEv = cev;
	intLevelConstraint=levcon;
};

int SR::EventMatrix::TotalEventsPending() {
#ifdef _DEBUG
	int rtnint=0;
	for (int i=0;i<vecSizeOfStacks.size();++i) rtnint+=vecSizeOfStacks[i];
	if (rtnint != intTotalEvents) SR::srerror("Prob.");
	return rtnint;
#else
	return intTotalEvents;
#endif
};

int SR::EventMatrix::EventsPendingAtThisTimeStep() {
	return vecSizeOfStacks[intCurrentMatrixTime];
};

void SR::EventMatrix::ClearAllEvents() {
	for (int i=0;i<intMaxDelay;++i) {
		vecSizeOfStacks[i]=0;
                vecNumberOfConstrainedEvents[i]=0;
	}
	intTotalEvents=0;
};

void SR::EventMatrix::ApplyEvents(SR::ParameterSet& p) {
	SR::EventMatrix* ptEm = this;
	// int debug=0;
	bool eventResult;
	vector<SR::UntimedEvent>::iterator ptEv = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents;
	vector<SR::UntimedEvent>::iterator ptEvEnd = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime];
	while (ptEv != ptEvEnd) {
		eventResult = ptEv->ue(ptEv->ptNode1,ptEv->ptNode2,*ptEm,p);
		ptEvEnd = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime];
		ptEv++;
		intTotalEvents--;
		// debug++;
	}
	vecSizeOfStacks[intCurrentMatrixTime]=0;
	vecNumberOfConstrainedEvents[intCurrentMatrixTime]=0;
	intCurrentMatrixTime++;
	if (intCurrentMatrixTime==intMaxDelay) intCurrentMatrixTime=0;
};

void SR::EventMatrix::ApplyEvents(SR::ParameterSet& p, SR::GatherRunInformation& ri) {
	static int charbefore;
	SR::EventMatrix* ptEm = this;
	bool eventResult;
	vector<SR::UntimedEvent>::iterator ptEv = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents;
	vector<SR::UntimedEvent>::iterator ptEvEnd = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime];
	while (ptEv != ptEvEnd) {
		// cerr << ptEv << "\t" << ptEv->ue << "\t" << ptEv->ptNode1 << "\t" << ptEv->ptNode2 << "\t" << ptEv-matStacksLinear.begin() << "\t" << matStacksLinear.end()-ptEv << "\n";		charbefore = ptEv->ptNode2->GetCharacteristic();
		eventResult = ptEv->ue(ptEv->ptNode1,ptEv->ptNode2,*ptEm,p);
		if (eventResult) {
			ri.RegisterEventAfterApplication(ptEv,charbefore);
		}
		ptEvEnd = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime];
		intTotalEvents--;
		ptEv++;
	}
	vecSizeOfStacks[intCurrentMatrixTime]=0;
	vecNumberOfConstrainedEvents[intCurrentMatrixTime]=0;
	intCurrentMatrixTime++;
	if (intCurrentMatrixTime==intMaxDelay) intCurrentMatrixTime=0;
};

string SR::EventMatrix::StringOutputTheStack() {
	ostringstream oss;
	for (int i=0;i<intMaxDelay;++i) {
		for (int j=0;j<vecSizeOfStacks[i];++j) {
			oss << "Stack\t" << matStacksLinear[i*intMaxEvents+j].ptNode1 << "\t" << matStacksLinear[i*intMaxEvents+j].ptNode2 << "\t" << matStacksLinear[i*intMaxEvents+j].ue << endl;
		}
	}
	return oss.str();
};

void SR::EventMatrix::RemoveDuplicatesAtCurrentTime() {
	// Need to implement this.
	// Otherwise parameterisation will be way off for workplace networks.
};

bool SR::EventMatrix::CurrentTimeHasEvents() {
	if (vecSizeOfStacks[intCurrentMatrixTime]==0) return false;
	return true;
};

SR::UntimedEvent SR::EventMatrix::PopTop() {
	if (vecSizeOfStacks[intCurrentMatrixTime]==0) SR::srerror("PopTop tried on null stack.\n");
	vecSizeOfStacks[intCurrentMatrixTime]--;
	intTotalEvents--;
	return matStacksLinear[intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime]];
};

void SR::EventMatrix::ExportEvents(SR::EventMatrix& targetEm) {
	vector<SR::UntimedEvent>::iterator ptEv,ptEvEnd;
	for (int i=0;i<intMaxDelay;++i) {
		ptEv = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents;
		ptEvEnd = matStacksLinear.begin()+intCurrentMatrixTime*intMaxEvents+vecSizeOfStacks[intCurrentMatrixTime];
		while (ptEv != ptEvEnd) {
			targetEm.AddEvent(*ptEv,i);
			intTotalEvents--;
			ptEv++;
		}
		vecSizeOfStacks[intCurrentMatrixTime]=0;
		vecNumberOfConstrainedEvents[intCurrentMatrixTime]=0;
		intCurrentMatrixTime++;
		if (intCurrentMatrixTime==intMaxDelay) intCurrentMatrixTime=0;
	}
};

void SR::EventMatrix::RemoveEventAtRandom(SR::ParameterSet& p) {
	static vector<SR::UntimedEvent>::iterator ptEvTop,ptEvDelete;
	static int indexRemoval,indexCurrent=0;
	if (intTotalEvents == 0) SR::srerror("void SR::EventMatrix::RemoveEventAtRandom(SR::ParameterSet& p)");
//	indexRemoval = static_cast<int>(static_cast<double>(intTotalEvents)*NR::ran1(p.intSeed));
	indexRemoval = static_cast<int>(static_cast<double>(intTotalEvents)*gsl_rng_uniform(glob_rng));
	while (indexRemoval > vecSizeOfStacks[indexCurrent]) {
		indexRemoval -= vecSizeOfStacks[indexCurrent];
		indexCurrent++;
	}
#ifdef _DEBUG
	if (vecSizeOfStacks[indexCurrent] == 0) SR::srerror("Ooops I did it again.");
#endif
	ptEvTop = matStacksLinear.begin()+indexCurrent*intMaxEvents+vecSizeOfStacks[indexCurrent]-1;
	ptEvDelete = matStacksLinear.begin()+indexCurrent*intMaxEvents+indexRemoval;
	swap(*ptEvTop,*ptEvDelete);
	vecSizeOfStacks[indexCurrent]--;
	intTotalEvents--;
};

void SR::EventMatrix::SaveStackToFile(string filename, SR::GridHex& gh, int (*f1)(SR::UNTIMEDEVENT ue)) {
	static int delay,intN1,intN2,intEv;
	ofstream ofs;
	ofs.open(filename.c_str());
	if (ofs.fail()) SR::srerror("Problem opening Event file for output");
	for (int i=0;i<intMaxDelay;++i) {
		if (i<intCurrentMatrixTime) delay = i+intMaxDelay-intCurrentMatrixTime;
		else delay = i-intCurrentMatrixTime;
		for (int j=0;j<vecSizeOfStacks[i];++j) {
			intN1 = matStacksLinear[i*intMaxEvents+j].ptNode1 - gh.FirstNode();
			intN2 = matStacksLinear[i*intMaxEvents+j].ptNode2 - gh.FirstNode();
			intEv = f1(matStacksLinear[i*intMaxEvents+j].ue);
			ofs << delay << "\t" << intN1 << "\t" << intN2 << "\t" << intEv << "\t" << "\n";
		}
	}
	ofs.close();
};


void SR::EventMatrix::LoadPendingEventsFromFile(string filename, SR::GridHex& gh, SR::UNTIMEDEVENT (*f1)(int)) {
	static int delay,intN1,intN2,intEv;
	// static SR::Node *ptN1, *ptN2;
	static SR::UntimedEvent ue;
	ifstream ifs;
	ifs.open(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem opening Event file for input");
	while (ifs >> delay >> intN1 >> intN2 >> intEv) {
		ue.ptNode1 = gh.FirstNode() + intN1;
		ue.ptNode2 = gh.FirstNode() + intN2;
		ue.ue = f1(intEv);
		AddEvent(ue,delay);
	}
	ifs.close();
};
