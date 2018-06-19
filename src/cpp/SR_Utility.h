#ifndef SR_INC_UTILITY
#define SR_INC_UTILITY

#include<iostream>
// #include"nr.h"
#include <fstream>
#include <complex>
#include"SR_Parameter.h"
#include<gsl/gsl_rng.h>

/*
[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {

	// Some useful stuff
	bool AppendStringToFile(string filename, string data);
	bool AppendStringToFileWithWaitIfRequired(string filename, string data);
	bool OpenNullFile(string filename,string headers);
	bool TestFilePresent(string filename);
	void AddToListExact(vector<pair<double, double> >& v, pair<double,double>& d, double epsilon);
	void AverageOfList(vector<pair<double, double> >& v, double d);
	string PrintList(vector<pair<double, double> >& v);
	string PrintListNoCR(vector<pair<double, double> >& v);
	bool pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
	template<class T> bool IsIn(vector<T>& vec, T val);
	template<class T> ofstream& operator<<(ofstream& ofs, vector<T>& v);
	template<class T> vector<T> ReadVectorAsBinaryFromFile(ifstream& ifs);
	template<class T> void BinWrite(ofstream& ofs, T v);
	template<class T> T BinRead(ifstream& ifs);
	void srerror(const string etext);
	inline double Distance(double x1, double y1, double x2, double y2, double xmax, double ymax) {

		// rectangular in kilometers giving km distance
		/*
		static double dx,dy;
		dx = x1-x2;dy=y1-y2;
		return sqrt(dx*dx+dy*dy);
		*/

		// polar in decimal long (x) lat (y) giving km distance
		static double dist;
		dist = 6378.7*acos(sin(y1)*sin(y2)+cos(y1)*cos(y2)* cos(x2-x1));
		if (dist > -1e100 && dist < 1e100) return dist;
		else return 0;

	};

	inline double Distance(double x1, double y1, double x2, double y2) {
		return Distance(x1,y1,x2,y2,1e100,1e100);
	};
	struct IntCoord {int x,y;};

	class Node;
	class EventMatrix;
	class ParameterSet;
	class GatherRunInformation;

	// Type defs
	typedef double (*ONE_D_DISTRIBUTION)(double mean, SR::ParameterSet &p);
	typedef double (*KERNEL)(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset);
	typedef bool (*UNTIMEDEVENT)(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p);
	typedef void (*PROCESS)(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em);
	struct UntimedEvent {SR::Node *ptNode1,*ptNode2;UNTIMEDEVENT ue;};

	bool eqStrictUntimedEvent(const SR::UntimedEvent& ue1, const SR::UntimedEvent& ue2);

}

template<class T> ofstream& SR::operator<<(ofstream& ofs, vector<T>& v) {
	// Only for classes where sizeof(T) is defined;
	static char *filePointer;
	static int size;
	size = v.size();
	filePointer = (char*)(&size);
	ofs.write(filePointer,sizeof(int));
	filePointer = (char*)(v.begin());
	ofs.write(filePointer,size*(sizeof(int)));
	return ofs;
};

template<class T> vector<T> SR::ReadVectorAsBinaryFromFile(ifstream& ifs) {
	static int size;
	static char *filePointer;
	filePointer = (char*)(&size);
	ifs.read(filePointer,sizeof(int));
	vector<T> rtnvec(size);
	filePointer = (char*)(rtnvec.begin());
	ifs.read(filePointer,size*sizeof(T));
	return rtnvec;
};

template<class T> bool SR::IsIn(vector<T>& vec, T val) {
	typename vector<T>::iterator pt = vec.begin();
	while (pt!=vec.end()) {
		if (*pt == val) return true;
	}
	return false;
};

template<class T> void SR::BinWrite(ofstream& ofs, T v) {
	static char *p;
	p = (char*)(&v);
	// cout << sizeof(v) << "\n";
	ofs.write(p,sizeof(v));
};

template<class T> T SR::BinRead(ifstream& ifs) {
	static char *p;
	T r;
	p = (char*)(&r);
	ifs.read(p,sizeof(r));
	return r;
};

#endif
