#include"SR_DensityField.h"

SR::DensityField::DensityField(double d, int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in) {
	nox=nox_in+1;
	noy=noy_in+1;
	vals = new double[nox*noy]; for (int i=0;i<nox*noy;++i) vals[i]=0;
	minx=minx_in;
	miny=miny_in;
	maxx=maxx_in;
	maxy=maxy_in;
	stepx = static_cast<double>((maxx-minx)/(1.0*nox-1.0));
	stepy = static_cast<double>((maxy-miny)/(1.0*noy-1.0));
	for (int i=0;i<nox;++i) {
		for (int j=0;j<noy;++j) {
			vals[i*noy+j]=d;
		}
	}
	maxval = d;
};

SR::DensityField::DensityField(double f(double,double), int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in) {
	nox=nox_in+1;
	noy=noy_in+1;
	vals = new double[nox*noy]; for (int i=0;i<nox*noy;++i) vals[i]=0;
	minx=minx_in;
	miny=miny_in;
	maxx=maxx_in;
	maxy=maxy_in;
	stepx = static_cast<double>((maxx-minx)/(1.0*nox-1.0));
	stepy = static_cast<double>((maxy-miny)/(1.0*noy-1.0));
	for (int i=0;i<nox;++i) {
		for (int j=0;j<noy;++j) {
			vals[i*noy+j]=f(1.0*i*stepx+minx,10*j*stepy+miny);
		}
	}
	CalcMaxVal();
};

SR::DensityField::DensityField(string filename, double stepx_in, double stepy_in, double minx_in, double maxx_in, double miny_in, double maxy_in) {
	static double dblRT=1e-10;
	ifstream ifs;
	double tmpx,tmpy,tmpd,total=0;
	int xindex,yindex;
	string junk;
	nox=static_cast<int>(((maxx_in-minx_in)/stepx_in+1))+1;
	noy=static_cast<int>(((maxy_in-miny_in)/stepy_in+1))+1;
	vals = new double[nox*noy]; for (int i=0;i<nox*noy;++i) vals[i]=0;
	minx=minx_in;
	miny=miny_in;
	stepx=stepx_in;
	stepy=stepy_in;
	maxx=minx+(nox-1)*stepx;
	maxy=miny+(noy-1)*stepy;
	double *ptVecDbl;
	ifs.open(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem opening density file");
	ifs >> junk >> junk >> junk;
	while (ifs >> tmpx >> tmpy >> tmpd) {
		if (tmpy < 0) break;
		if (fmod(tmpx-minx,stepx)>dblRT || fmod(tmpy-miny,stepy)>dblRT) SR::srerror("Coords in density field file do not align with matrix");
		xindex = static_cast<int>((tmpx-minx)/stepx);yindex = static_cast<int>((tmpy-miny)/stepy);
		if (xindex < nox && xindex >= 0 && yindex < noy && yindex >= 0) {
			if (vals[xindex*noy+yindex]!=0) cerr << "Warning : Multiple assignments while reading density field\n";
			vals[xindex*noy+yindex]=tmpd;
			total+=tmpd;
		}
	};
	ifs.close();
	if (tmpy < 0) {
		cerr << "Uniform density selected with a -ve y coord in DensityField::DensityField\n";
		ptVecDbl = vals;
		while (ptVecDbl!=vals+nox*noy) {
			*ptVecDbl = 1;
			ptVecDbl++;
		}
	}
	CalcMaxVal();
	cerr << "Density field constructed.  Sum of densities (0 for uniform): " << total << "\n";
};

SR::DensityField::DensityField(string filename) {

	static int nullData;
	ifstream ifs;
	double tmp,total=0; //This is never increased?
	string junk;
	vector<double>::iterator ptVecDbl;
	ifs.open(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem opening density file");

	ifs >> junk >> nox;
	ifs >> junk >> noy;
	ifs >> junk >> minx;
	ifs >> junk >> miny;
	ifs >> junk >> stepx;
	ifs >> junk >> nullData;
	minx = minx / 57.2957795131;
	miny = miny / 57.2957795131;
	stepx = stepx / 57.2957795131;
	stepy=stepx;

	vals = new double[nox*noy]; for (int i=0;i<nox*noy;++i) vals[i]=0;

	maxx=minx+(nox-1)*stepx;
	maxy=miny+(noy-1)*stepy;

	for (int i=noy-1;i>=0;--i) {
		for (int j=0;j<nox;++j) {
			ifs >> tmp;
			if (tmp == nullData) tmp = 0;
			//Add total += tmp; here to track total?
			vals[j*noy+(noy-i-1)]=tmp;
		}
	}

	ifs.close();
	CalcMaxVal();
	cerr << "Density field constructed.  Sum of densities (0 for uniform): " << total << "\n";

};

string SR::DensityField::Table() {
	ostringstream oss;
	oss << "0\t";
	for (int i=0;i<nox;++i) oss << minx+1.0*i*stepx << "\t";
	oss << "\n";
	for (int j=0;j<noy;++j) {
		oss << 1.0*j*stepy+miny << "\t";
		for (int i=0;i<nox;++i) oss << vals[i*noy+j] << "\t";
		oss << "\n";
	}
	return oss.str();
};

double SR::DensityField::Value(double x, double y) {
	double rtnval;
	static int xcoord,ycoord;
	// cerr << "here7" << endl;
	if (x < minx || x > maxx || y < miny || y > maxy) SR::srerror("coords out of range in DensityField::Value(double x, double y)");
	// cerr << "here8" << endl;
	xcoord = static_cast<int>((x-minx)/stepx);
	ycoord = static_cast<int>((y-miny)/stepy);
	// cerr << "here9" << endl;
	// Consider edit here
	// rtnval = (vals[xcoord*noy+ycoord]+vals[(xcoord+1)*noy+ycoord]+vals[(xcoord+1)*noy+ycoord+1]+vals[xcoord*noy+ycoord+1])/4;
	rtnval = vals[xcoord*noy+ycoord];
	return rtnval;
};

void SR::DensityField::CalcMaxVal() {
	double* pt = vals;
	double rtnval = -9999;
	while (pt != vals+nox*noy) {
		if (*pt > rtnval) rtnval = *pt;
		pt++;
	}
	if (rtnval==0) SR::srerror("Entire density field is zero.");
	maxval = rtnval;
};

void SR::DensityField::WriteAsciiGrid(string filename) {

	static int nullData;
	ofstream ofs;
	double tmp,total=0;
	string junk;
	vector<double>::iterator ptVecDbl;
	ofs.open(filename.c_str());
	if (ofs.fail()) SR::srerror("Problem opening density file to write");

	ofs << "NCOLS" << " " << nox << "\n";
	ofs << "NROWS" << " " << noy << "\n";
	ofs << "XLLCORNER" << " " << minx*57.2957795131 << "\n";
	ofs << "YLLCORNER" << " " << miny*57.2957795131 << "\n";
	ofs << "CELLSIZE" << " " << stepx*57.2957795131 << "\n";
	ofs << "NODATA_value -2147483647\n";

	for (int i=noy-1;i>=0;--i) {
		for (int j=0;j<nox;++j) {
			ofs << vals[j*noy+(noy-i-1)] << " ";
		}
		ofs << "\n";
	}

	ofs.close();
	cerr << "Density written to field " << filename << "\n";

};

