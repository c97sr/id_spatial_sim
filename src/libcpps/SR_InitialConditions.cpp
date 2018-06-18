#include"SR_InitialConditions.h"

extern gsl_rng * glob_rng;

SR::InitialConditions::InitialConditions(GridHex& g, ParameterSet& p, int no, int msn, double x_in, double y_in, double r_in, int ms) :
vecNodes(msn,-1) {

	// static double r_dash;

	maxsamples=ms;
	x=x_in / 57.29577951,y=y_in / 57.29577951,r=r_in;
	intHexagonsUsed=0;
	intTotalNodes=0;
	currentsamples=0;
	initialnumber=no;

	// r_dash = r / 6378.7;

	//  Max no hex = 2 * pi / 3 / sqrt(3) * (R+w)^2 / w^2
	// Consider upping to (R+2w) if required
	double factor = 1.209199576;
	double hexwidth = g.GetHexagonWidth();
	int maxNoHexagons = static_cast<int>(factor*(r+2*hexwidth)*(r+2*hexwidth)/hexwidth/hexwidth);

	if (maxNoHexagons > g.GetNoHexagons()*0.5) blLocalSeeding = false;
	else blLocalSeeding = true;

	SR::Hexagon* ptCurrentHexagon = g.FirstHexagon();

	if (blLocalSeeding) {

		ptvecListOfHexagons = new int[maxNoHexagons];
		ptvecNodesInHexagons = new int[maxNoHexagons];

		//  Loop to add in all hexagons with centres inside the r+w width
		while (ptCurrentHexagon != g.LastHexagon()) {
			if (SR::Distance(ptCurrentHexagon->GetX(),ptCurrentHexagon->GetY(),x,y,1e10,1e10) <= r + hexwidth) {
				ptvecListOfHexagons[intHexagonsUsed] = ptCurrentHexagon - g.FirstHexagon();
				ptvecNodesInHexagons[intHexagonsUsed] = ptCurrentHexagon->GetNoNodes();
				intHexagonsUsed++;
				intTotalNodes+=ptCurrentHexagon->GetNoNodes();
				if (intHexagonsUsed == maxNoHexagons)
					SR::srerror("Maximum number of seeding hexagons incorrectly calculated.");
			}
			ptCurrentHexagon++;
		}
	}
};

SR::InitialConditions::~InitialConditions() {
	if (blLocalSeeding) {
		delete [] ptvecListOfHexagons;
		delete [] ptvecNodesInHexagons;
	}
};

void SR::InitialConditions::Reselect(GridHex& g, ParameterSet& p) {
	// static int number_requested;
	Reselect(g,p,initialnumber);
};

void SR::InitialConditions::Reselect(GridHex& g, ParameterSet& p, int number) {
	double nonodes = static_cast<double>(g.GetNoNodes());
	int currenttrys=0,currenthits=0;
	int intCurrentHex,intSelectedNode;
	bool notused;
	if (number + currentsamples > static_cast<int>(vecNodes.size())) SR::srerror("Too many seeds requested. Is max seed number too high?");
	SR::Node *ptNode;
	while (currenthits < number && currenttrys < maxsamples) {
		if (blLocalSeeding) {
			intCurrentHex=0;
//			intSelectedNode = static_cast<int>(NR::ran2(p.intSeed)*intTotalNodes);
			intSelectedNode = static_cast<int>(gsl_rng_uniform(glob_rng)*intTotalNodes);
			while (intSelectedNode >= ptvecNodesInHexagons[intCurrentHex]) {
				intSelectedNode-=ptvecNodesInHexagons[intCurrentHex];
				intCurrentHex++;
				if (intCurrentHex==intHexagonsUsed)
					SR::srerror("Problem with local seed selection.");
			}
			ptNode = *((g.FirstHexagon()+ptvecListOfHexagons[intCurrentHex])->GetFirstNode()+intSelectedNode);
		} else {
//			ptNode = g.FirstNode()+static_cast<int>(NR::ran2(p.intSeed)*nonodes);
			ptNode = g.FirstNode()+static_cast<int>(gsl_rng_uniform(glob_rng)*nonodes);
		}
		if (SR::Distance(ptNode->GetX(),ptNode->GetY(),x,y,1e10,1e10) <= r)  {
			vector<int>::iterator ptInt = vecNodes.begin();
			notused = true;
			while (ptInt != vecNodes.begin()+currenthits && notused) {
				if (ptNode->GetIndex()==*ptInt) notused=false;
				ptInt++;
			}
			if (notused) {
				*ptInt=ptNode->GetIndex();
				currenthits++;
				currentsamples++;
			}
		}
		currenttrys++;
	};
	if (currenttrys==maxsamples)
		SR::srerror("Too many samples taken in setting up or reselecting initial conditions.");
};

void SR::InitialConditions::ApplySeed(GridHex& g, SR::ParameterSet& p, PROCESS proc, EventMatrix& em) {
	vector<int>::iterator itInt = vecNodes.begin();
	while (currentsamples > 0) {
		proc(p,g.LastNode()-1,g.FirstNode()+*itInt,em);
		itInt++;
		currentsamples--;
	}
};

void SR::InitialConditions::Trickle(GridHex& g, SR::ParameterSet& p, PROCESS proc, EventMatrix& em, double trate, double tdur) {
	static int intNumberRequired;
	static double *t;
	t = p.GetPointer("dblCurrentTime");
//	if (trate > 0 && tdur > *t) intNumberRequired = static_cast<int>(NR::poidev(trate,p.intSeed));
	if (trate > 0 && tdur > *t) intNumberRequired = static_cast<int>(gsl_ran_poisson(glob_rng,trate));
	else intNumberRequired=0;
	if (intNumberRequired > maxsamples) intNumberRequired = maxsamples;
	Reselect(g,p,intNumberRequired);
	ApplySeed(g,p,proc,em);
};

