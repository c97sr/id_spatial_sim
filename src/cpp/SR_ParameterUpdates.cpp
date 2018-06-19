#include"SR_ParameterUpdates.h"

extern gsl_rng * glob_rng;

SR::ParameterUpdates::ParameterUpdates(SR::ParameterSet& ps_in, string inFileName, string outFileName, string uniqueId_in) :
paramNames(mp), paramPointers(mp), paramInitial(mp), paramCurrent(mp), paramWeight(mp),
paramMin(mp), paramMax(mp), paramPropStep(mp),   paramRangeIsLog(mp) {
	noParams=0;
	sampleNo=0;
	currentLnLike = -1e100;
	dblLogProposalFactor = -1e100;
	intNoChosenParam=0;
	totalWeight=0;
	ps=&ps_in;
	// Parameter file needs to be structured in correct way
	// [name		initial	min	max		weight	step	log]
	string tmp, filecontents;
	ifstream ifs;
	ostringstream oss;
	outputFile.Set(outFileName);
	inputFile.Set(inFileName);
	uniqueId.Set(uniqueId_in);
	ifs.open(inputFile.Get().c_str());
	if (ifs.fail()) SR::srerror("Problem with the input parameter file.");
	while (ifs >> tmp) filecontents += tmp + "\t";
	ifs.close();
	double dblTmp;
	bool blTmp;
	filecontents = SR::StripComments(filecontents);
	istringstream iss(filecontents);
	// cerr << "\n" << filecontents << "\n";
	while (iss >> tmp) {
		if (noParams == mp) SR::srerror("Too many fitted parameters in SR::ParameterUpdates::ParameterUpdates");
		(paramNames.begin()+noParams)->Set(tmp);
		iss >> dblTmp; paramInitial[noParams]=dblTmp;
		iss >> dblTmp; paramMin[noParams]=dblTmp;
		iss >> dblTmp; paramMax[noParams]=dblTmp;
		iss >> dblTmp; paramWeight[noParams]=dblTmp;
		iss >> dblTmp; paramPropStep[noParams]=dblTmp;
		iss >> blTmp; paramRangeIsLog[noParams]=blTmp;
		noParams++;
	}
	if (noParams==0) SR::srerror("Null sampling file in ParameterUpdates constructor.");
	// Set up the pointers
	for (int i=0;i<noParams;++i) {
		paramPointers[i] = ps->GetPointer(paramNames[i].Get());
		*(paramPointers[i]) = paramInitial[i];
		totalWeight+=paramWeight[i];
	}
	// Validate output file
	if (!TestFilePresent(outputFile.Get())) {
		cerr << "Output file not present for parameter updates.  Initialising: " << outFileName << "...";
		oss << "";
		for (int i=0;i<noParams;++i) {
			oss << paramNames[i].Get() << "\t" ;
		}
		oss << "lnlike";
		if (!SR::OpenNullFile(outFileName,oss.str()))
			SR::srerror("Problem opening file in parameter update constructor.");
		cerr << "done.\n";
	}
};

void SR::ParameterUpdates::LogParameterValues(double lnlike, int count) {
	ostringstream oss;
	// oss << count << "\t";
	// for (int i=0;i<noParams;++i) oss << paramNames[i].Get() << ":\t" << *(paramPointers[i]) << "\t";
	for (int i=0;i<noParams;++i) oss << *(paramPointers[i]) << "\t";
	oss << lnlike << "\n";
	SR::AppendStringToFileWithWaitIfRequired(outputFile.Get(),oss.str());
};

void SR::ParameterUpdates::LogParameterValues(double lnlike, ostringstream &oss, int count) {
	// oss << count << "\t";
	// for (int i=0;i<noParams;++i) oss << paramNames[i].Get() << ":\t" << *(paramPointers[i]) << "\t";
	for (int i=0;i<noParams;++i) oss << *(paramPointers[i]) << "\t";
	oss << lnlike << "\n";
};

void SR::ParameterUpdates::ProposeUpdate() {
	static double currentWeight,randomWeight,dblParamValue,dblJump,dblMaxJump;
	static double *param;

	int dbgtest = r->val();

//	randomWeight = totalWeight*NR::ran2(ps->intSeed);
	randomWeight = totalWeight*gsl_rng_uniform(glob_rng);
	intNoChosenParam = 0;
	currentWeight = paramWeight[0];
	while (randomWeight > currentWeight) {
		intNoChosenParam++;
		currentWeight += paramWeight[intNoChosenParam];
	}

	// Transform onto the 0-1 uniform scale
	param = paramPointers[intNoChosenParam];
	paramCurrent[intNoChosenParam]=*param;

	if (paramRangeIsLog[intNoChosenParam])
		dblParamValue = 1.0 / (log(paramMax[intNoChosenParam])-log(paramMin[intNoChosenParam])) *
		log((*param) / paramMin[intNoChosenParam]);
	else dblParamValue = (*param - paramMin[intNoChosenParam]) /
		(paramMax[intNoChosenParam] - paramMin[intNoChosenParam]);

	// Calculate the new value for the variable
	dblMaxJump = paramPropStep[intNoChosenParam];
//	dblJump = NR::ran2(ps->intSeed)*dblMaxJump-dblMaxJump/2.0;
	dblJump = gsl_rng_uniform(glob_rng)*dblMaxJump-dblMaxJump/2.0;
	dblParamValue += dblJump;
	if (dblParamValue < 0) dblParamValue= -1.0*dblParamValue;
	if (dblParamValue >= 1) dblParamValue=2-dblParamValue;

	// Transform back to the appropriate scale
	if (paramRangeIsLog[intNoChosenParam]) {
		*param = paramMin[intNoChosenParam]*pow(10.0,(log10(paramMax[intNoChosenParam])-log10(paramMin[intNoChosenParam]))*dblParamValue);
	} else {
		*param = paramMin[intNoChosenParam]+(paramMax[intNoChosenParam]-paramMin[intNoChosenParam])*dblParamValue;
	}
};

void SR::ParameterUpdates::AcceptUpdate() {
	intNoChosenParam = -1;
};

void SR::ParameterUpdates::DeclineUpdate() {
	*(paramPointers[intNoChosenParam])=paramCurrent[intNoChosenParam];
	intNoChosenParam = -1;
};

