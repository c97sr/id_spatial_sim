#include"SR_Utility.h"

bool SR::eqStrictUntimedEvent(const SR::UntimedEvent& ue1, const SR::UntimedEvent& ue2) {
	if (ue1.ptNode1==ue2.ptNode1 && ue1.ptNode2==ue2.ptNode2 && ue1.ue == ue2.ue) return true;
	else return false;
};

bool SR::AppendStringToFile(string filename, string data) {
	ofstream s(filename.c_str(),iostream::app);
	s.precision(6);
	if (s.fail()) return false;
	s << data;
	s.close();
	return true;
};

bool SR::AppendStringToFileWithWaitIfRequired(string filename, string data) {
	static int seed = 1234;
	static int maxtries = 10;
	int currentfails = 0;
	double maxwaittime=10*CLOCKS_PER_SEC;
	double waittime;
	clock_t delay;
	while (!SR::AppendStringToFile(filename,data) && currentfails != maxtries) {
		waittime = NR::ran2(seed)*maxwaittime;
		cerr << "Failed on attempt " << currentfails+1 <<  " to append to " << filename;
		delay = static_cast<long unsigned int>(waittime)+clock();
		while (delay > clock());
		currentfails++;
	};
	if (currentfails==maxtries) return false;
	else return true;
};

bool SR::OpenNullFile(string filename, string headers) {
	ofstream s(filename.c_str());
	if (s.fail()) return false;
	s << headers << "\n";
	s.close();
	return true;
};

bool SR::TestFilePresent(string filename) {
	ifstream s;
	s.open(filename.c_str());
	bool rtn;
	if (s.fail()) rtn = false;
	else rtn = true;
	s.close();
	return rtn;
};


void SR::AddToListExact(vector<pair<double, double> >& v, pair<double,double>& d, double epsilon) {
	vector<pair<double, double> >::iterator it = v.begin();
	while (it != v.end()) {
		if (fabs(it->first - d.first) < epsilon) {
			it->second += d.second;
			return;
		}
		it++;
	}
	SR::srerror("Value not found in void AddToListExact(vector<pair<double, double> >& v, pair<double, double>& d) ");
}

void SR::AverageOfList(vector<pair<double, double> >& v, double d) {
	vector<pair<double, double> >::iterator it = v.begin();
	while (it != v.end()) {
		it->second /= d;
		it++;
	}
}

string SR::PrintList(vector<pair<double, double> >& v) {
	ostringstream oss;
	vector<pair<double, double> >::iterator it = v.begin();
	while (it != v.end()) {
		oss << it->first << "\t" << it->second << "\n";
		it++;
	}
	return oss.str();
}

string SR::PrintListNoCR(vector<pair<double, double> >& v) {
	ostringstream oss;
	vector<pair<double, double> >::iterator it = v.begin();
	while (it != v.end()) {
		oss << it->second << "\t";
		it++;
	}
	return oss.str();
}

bool SR::pnpoly(int nvert, double *vertx, double *verty, double testx, double testy) {
	// Adapted from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
	int i, j;
	bool c = false;
	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((verty[i]>testy) != (verty[j]>testy)) &&
				(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
			c = !c;
	}
  return c;
}

void SR::srerror(const string etext) {
	cerr << "Run-time error..." << endl;
	cerr << etext << endl;
	cerr << "...now exiting." << endl;
	exit(1);
};



