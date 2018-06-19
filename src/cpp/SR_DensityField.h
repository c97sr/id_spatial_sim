#ifndef SR_INC_DENSITYFIELD
#define SR_INC_DENSITYFIELD

#include<iostream>
#include<vector>
// #include"nr.h"
#include"SR_Utility.h"

/*

[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {

	// inline void srerror(const string etext);

	class DensityField {
	private:
		double* vals;
		double minx,miny,maxx,maxy,stepx,stepy,maxval;
		int nox,noy; // number of grid points
		void CalcMaxVal();
	public:
		// nox and noy are numbers of intervals in the constructors, not number of points
		DensityField() {SR::srerror("No default constructor for DensityField");};
		DensityField(double d, int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
		DensityField(double f(double,double), int nox_in, int noy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
		DensityField(string filename, double stepx_in, double stepy_in, double minx_in, double maxx_in, double miny_in, double maxy_in);
		DensityField(string filename);
		~DensityField(){delete [] vals;};
		double Value(double x, double y);
		string Table();
		inline double GetMinX(){return minx;};
		inline double GetMaxX(){return maxx;};
		inline double GetMinY(){return miny;};
		inline double GetMaxY(){return maxy;};
		inline double GetMaxVal(){return maxval;};
		void WriteAsciiGrid(string filename);
	};
}

#endif
