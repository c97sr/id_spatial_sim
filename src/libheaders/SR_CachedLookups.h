#ifndef SR_INC_CACHEDLOOKUPS
#define SR_INC_CACHEDLOOKUPS

#include<iostream>
#include<vector>
#include"nr.h"
#include"SR_Parameter.h"

/*
[Enter details]
S Riley [Enter details]
*/

using namespace std;

namespace SR {

	// inline void srerror(const string etext);

	class CachedIntLookup {
	private:
		double* values;
		int intNoValues;
		double dx;
		int min,max;
		double* ptStart;
	public:
		CachedIntLookup() {SR::srerror("No default constructor for CachedIntLookup.");};
		CachedIntLookup(double f(int,SR::ParameterSet&),int min_in,int max_in, SR::ParameterSet& p);
		~CachedIntLookup();
		double GetCheckedValue(int i);
		double GetFastValue(int i);
	};

	class CachedDblLookup {
	private:
		static const int max_values = 10000;
		int values_used;
		double* values;
		double dlnx,start,logstart,max;
		double* ptStart;
		bool BeenMadeHistogram;
	public:
		CachedDblLookup() {SR::srerror("No default constructor for CachedIntLookup.");};
		CachedDblLookup(double f(double,SR::ParameterSet&),double start, double incs_per_order, double no_orders, SR::ParameterSet& p);
		CachedDblLookup(double f(double,int, double*),double start_in, double incs_per_order, double no_orders,int np, double *vecp);
		CachedDblLookup(string filename);
		~CachedDblLookup();
		double GetCheckedValue(double d);
		double GetFastValue(double d);
		void Normalise();
		void MakeHistogram();
		double GetHistValue(double d);
		string PrintDistributionOfCommutes();
	};
};

#endif
