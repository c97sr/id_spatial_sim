#include"SR_Parameter.h"
#include"SR_Utility.h"

void SR::ParameterDash::Assign(string s, bool n, double d, string t) {
	name.Set(s);
	numeric=n;
	value=d;
	tag.Set(t);
}

vector<SR::ParameterDash>::iterator SR::ParameterSet::FindName(string s){
	vector<SR::ParameterDash>::iterator it = p.begin();
	vector<SR::ParameterDash>::iterator e = it + n;
	while (it!=e && it!=p.end()) {
		if (it->GetName()==s) return it;
		++it;
	}
	if (it==p.end()) SR::srerror("Not enough space for a new parameter name.");
	return it;
};

void SR::ParameterSet::AddValue(string s, double d) {
	if (locked) 
		SR::srerror("Not allowed to add value if parameters are locked.");
	vector<ParameterDash>::iterator it = FindName(s);
	double tmp = it->Value();
	it->Set(d);
	if (it-p.begin()<n) {
		if (!it->IsNumeric()) SR::srerror("Expecting a numeric value in void SR::ParameterSet::AddValue(string s, double d)");
		cerr << "Overwriting value of " << s << ". Old value: " << tmp << ". New value: " << d << "\n";
	} else {
		it->SetNumeric(true);
		it->SetName(s);
		n++;
	};
};

void SR::ParameterSet::AddValue(string s, string t) {
	if (locked) SR::srerror("Not allowed to add value if parameters are locked.");
	vector<ParameterDash>::iterator it = FindName(s);
	it->Set(t);
	if (it-p.begin()<n) {
		if (it->IsNumeric()) SR::srerror("Expecting an alpha value in void SR::ParameterSet::AddValue(string s, double d)");
		cerr << "Overwriting value of " << s << ". Old value: " << it->Value() << ". New value: " << t << "\n";
	} else {
		it->SetNumeric(false);
		it->SetName(s);
		n++;
	};
};

double* SR::ParameterSet::GetPointer(string s) {
	vector<ParameterDash>::iterator it = FindName(s);
	if (!it->IsNumeric()) SR::srerror("Pointers only to numerics in SR::ParameterSet::GetPointer(string s)");
	if (it-p.begin()>=n) SR::srerror(("Problem in SR::ParameterSet::GetPointer(string s): "+s).c_str());
	return it->Pointer();
};

double* SR::ParameterSet::GetInitialPointer(string s) {
	vector<ParameterDash>::iterator it = FindName(s);
	if (!it->IsNumeric()) SR::srerror("Pointers only to numerics in SR::ParameterSet::GetPointer(string s)");
	if (it-p.begin()>=n) SR::srerror(("Problem in SR::ParameterSet::GetPointer(string s): "+s).c_str());
	return it->InitialPointer();
};

vector<SR::ParameterDash>::iterator SR::ParameterSet::GetParameterPointer(string s) {
	vector<ParameterDash>::iterator it = FindName(s);
	if (!it->IsNumeric()) SR::srerror("Pointers only to numerics in SR::ParameterSet::GetPointer(string s)");
	if (it-p.begin()>=n) SR::srerror(("Problem in SR::ParameterSet::GetPointer(string s): "+s).c_str());
	return it;
};

double SR::ParameterSet::GetValue(string s) {
	vector<ParameterDash>::iterator it = FindName(s);
	if (it-p.begin()>=n) SR::srerror(("Problem in SR::ParameterSet::GetValue(string s): "+s).c_str());
	if (!it->IsNumeric()) SR::srerror("Value not defined in SR::ParameterSet::GetPointer(string s)");
	return it->Value();
};

string SR::ParameterSet::GetTag(string s) {
	vector<ParameterDash>::iterator it = FindName(s);
	if (it-p.begin()>=n) SR::srerror(("Problem in SR::ParameterSet::GetValue(string s): "+s).c_str());
	if (it->IsNumeric()) SR::srerror("Tag not defined in SR::ParameterSet::GetTag(string s)");
	return it->Tag();
};

void SR::ParameterSet::ReadParams(string s) {
	istringstream iss(s);
	string name,tag;
	bool n,seednotpresent=true;
	double val;
	while (iss >> name) {
		iss >> n;
		if (name=="intSeed") {
			if (!n) SR::srerror("intSeed must be an int.");
			iss >> intSeed;
			if (intSeed > -1) intSeed *= -1;
			seednotpresent=false;
		} else {
			if (n) {
				iss >> val;
				AddValue(name,val);
			} else {
				iss >> tag;
				AddValue(name,tag);
			}
		}
	}
	// Every parameter file needs to have an intSeed
	// if (seednotpresent) cerr << "Warning: seed not initialised during ReadParams\n";	
};

string SR::ParameterSet::WriteParams() {
	vector<ParameterDash>::iterator it;
	ostringstream oss;
	it = p.begin();
	while (it != p.begin()+n) {
		oss << it->GetName() << "\t" << it->IsNumeric() << "\t";
		if (it->IsNumeric()) oss << it->Value();
		else oss << it->Tag();
		oss << "\n";
		it++;
	}
	return oss.str();
};

void SR::ParameterSet::Lock() {
	vector<SR::ParameterDash>::iterator it = p.begin();
	while (it != p.begin()+n) {
		it->SetInitialValue();
		it++;
	}
	locked=true;
};

void SR::ParameterSet::RevertToInitialValues(){
	vector<SR::ParameterDash>::iterator it = p.begin();
	while (it != p.begin()+n) {
		it->RevertToInitialValue();
		it++;
	}
	locked=true;
};

ofstream& SR::operator<<(ofstream& ofs, ParameterDash& p) {
	static char *filePointer;
	filePointer = (char*)(&p.numeric); ofs.write(filePointer,sizeof(bool));
	ofs << p.name; 
	ofs << p.tag;
	filePointer = (char*)(&p.value); ofs.write(filePointer,sizeof(double));
	filePointer = (char*)(&p.initialvalue); ofs.write(filePointer,sizeof(double));
	return ofs;
};

ifstream& SR::operator>>(ifstream& ifs, ParameterDash& p) {
	static char *filePointer;
	filePointer = (char*)(&p.numeric); ifs.read(filePointer,sizeof(bool));
	ifs >> p.name;
	ifs >> p.tag;
	filePointer = (char*)(&p.value); ifs.read(filePointer,sizeof(double));
	filePointer = (char*)(&p.initialvalue); ifs.read(filePointer,sizeof(double));
	return ifs;
};

ofstream& SR::operator<<(ofstream& ofs, ParameterSet& ps) {
	static char *filePointer;
	filePointer = (char*)(&ps.n); ofs.write(filePointer,sizeof(int));
	filePointer = (char*)(&ps.locked); ofs.write(filePointer,sizeof(bool));
	for (int i=0;i<ps.max;++i) ofs << ps.p[i];
	return ofs;
};

ifstream& SR::operator>>(ifstream& ifs, ParameterSet& ps) {
	static char *filePointer;
	filePointer = (char*)(&ps.n); ifs.read(filePointer,sizeof(int));
	filePointer = (char*)(&ps.locked); ifs.read(filePointer,sizeof(bool));
	for (int i=0;i<ps.max;++i) ifs >> ps.p[i];
	ps.locked = false;
	return ifs;
};	

void SR::ParameterSet::ReadParamsFromFile(string s) {
	string filecontents,tmp;
	ifstream ifs(s.c_str());
	if (ifs.fail()) SR::srerror("Problem with parameter file.");
	while (ifs >> tmp) filecontents += tmp + "\t";
	ifs.close();
	filecontents=SR::StripComments(filecontents);
	ReadParams(filecontents);
};

SR::ParameterValueSet::ParameterValueSet(string filename) {
	int norows=0;
	string tmp;
	double dbltmp;
	int inttmp;
	maxsetnumber=0;
	ifstream ifs(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem with parameter file:SR::ParameterValueSet::ParameterValueSet(string filename)");
	while (ifs >> tmp) {
		ifs >> tmp; ifs >> tmp; 
		norows++;
	}
	ifs.close();

	noChanges = norows;
	norows=0;

	int dbgtmp = noChanges;

	paramlabels = new (string(*[dbgtmp]));
	paramvalues = new double[noChanges];
	setnumbers = new int[noChanges];

	ifstream ifs2(filename.c_str());
	if (ifs2.fail()) SR::srerror("Problem with parameter file:SR::ParameterValueSet::ParameterValueSet(string filename). 2nd call");
	while (ifs2 >> inttmp) {
		setnumbers[norows]=inttmp;
		if (inttmp > maxsetnumber) maxsetnumber=inttmp;
		ifs2 >> tmp; *(paramlabels+norows) = new string; **(paramlabels+norows)=tmp; 
		ifs2 >> dbltmp; paramvalues[norows]=dbltmp;
		norows++;
	}	
	if (norows != noChanges) SR::srerror("Something wrong in SR::ParameterValueSet::ParameterValueSet(string filename)");
	ifs2.close();
};

SR::ParameterValueSetB::ParameterValueSetB(string filename) {
	string tmp;
	bool loopNotStop;
	double dbltmp;

	noParams=0;
	noChanges=0;
	ifstream ifs(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem with parameter file:SR::ParameterValueSetB::ParameterValueSet(string filename)");
	loopNotStop=true;
	while (loopNotStop) {
		ifs >> tmp;
		if (tmp=="PARAM_LINE_END") loopNotStop = false;
		else noParams++;
	}
	loopNotStop=true;
	while (ifs >> tmp) {
		for (int i=0;i<noParams-1;++i) ifs >> tmp;
		noChanges++;
	}
	ifs.close();

	int dbgtmp = noParams;

	paramlabels = new (string(*[dbgtmp]));
	paramvalues = new double[noChanges*noParams];

	ifstream ifs2(filename.c_str());
	if (ifs2.fail()) SR::srerror("Problem with parameter file:SR::ParameterValueSet::ParameterValueSet(string filename). 2nd call");
	for (int i=0;i<noParams;++i) {ifs2 >> tmp;*(paramlabels+i) = new string; **(paramlabels+i)=tmp;}
	ifs2 >> tmp;
	for (int i=0;i<noChanges;++i) {
		for (int j=0;j<noParams;++j) {ifs2 >> dbltmp;paramvalues[i*noParams+j]=dbltmp;}
	}
	ifs2.close();
	currentparamnumber=0;
};

SR::ParameterValueSet::~ParameterValueSet() {
	for (int i=0;i<noChanges;++i) delete *(paramlabels+i); 
	delete [] paramlabels;
	delete [] paramvalues;
	delete [] setnumbers;
};

SR::ParameterValueSetB::~ParameterValueSetB() {
	for (int i=0;i<noParams;++i) delete *(paramlabels+i); 
	delete [] paramlabels;
	delete [] paramvalues;
};

void SR::ParameterValueSet::UpdateParameterSet(SR::ParameterSet &p, int ps_num) {
	static int *ptSetNo;
	static string **ptLabel;
	static double *val;
	ptSetNo = setnumbers;
	ptLabel = paramlabels;
	val = paramvalues;
	for (int i=0;i<noChanges;++i) {
		// Change line below to do actual change
		// if (*ptSetNo == ps_num)	cerr << "Trying to change " << **ptLabel << " to " << *val << "\n";
		if (*ptSetNo == ps_num)	p.ChangeValue(**ptLabel,*val);
		ptSetNo++;
		ptLabel++;
		val++;		
	}
}

void SR::ParameterValueSetB::UpdateParameterSet(SR::ParameterSet &p) {
	if (currentparamnumber==noChanges) SR::srerror("Attempted to overrun param changes.");
	for (int i=0;i<noParams;++i) {
		p.ChangeValue(**(paramlabels+i),paramvalues[currentparamnumber*noParams+i]);
	}
	currentparamnumber++;
}

SR::ParameterSet SR::CommandLineAndParameterFile(int argc, char **argv) {

	// Read from command line, put parameter file name into strParamFile and argument list into strArgs
	int intNoArgs = 1;
	if (argc<intNoArgs+1) SR::srerror("First argument must be name of parameter file. Rest parsed as parameter values.\n");
	string strParamFile, strArgs;
	strParamFile = argv[1];
	if ((argc-(intNoArgs+1))%3!=0) SR::srerror("An even number of parameter arguments are required.\n");
	for (int i=intNoArgs+1;i<argc;i=i+3){strArgs+=argv[i];strArgs+="\t";strArgs+=argv[i+1];strArgs+="\t";strArgs+=argv[i+2];strArgs+="\t";};

	// Read parameters from param file and then update with command line params
	SR::ParameterSet rtnval;
	rtnval.ReadParamsFromFile(strParamFile);
	rtnval.ReadParams(strArgs);

	return rtnval;
};

