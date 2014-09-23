#ifndef _SR_INC_HAZARD
#define _SR_INC_HAZARD

/*
XXXX

S Riley XXXX

*/

#include<iostream>
#include"nr.h"
#include"SR_Utility.h"
#include"SR_Parameter.h"

using namespace std;

namespace SR {
	class Hazard {
	private:
		int intNoPoints;
		double *ptFirstHazard;
		double *ptFirstLnHazard;
		double *ptFirstIntegratedHazard;
		double dblMin, dblMax, dblDT;
		double (*h)(SR::ParameterSet& p, double t);
		double LinearInterpolate(double *ptM, double t);
	public:
		Hazard(){SR::srerror("No default constructor for Hazard");};
		Hazard(double (*h_in)(SR::ParameterSet&, double), int resolution, double min_in, double max_in, SR::ParameterSet& p);
		~Hazard();
		void RecalculateHazard(SR::ParameterSet& p);
		void RecalculateLnHazard(SR::ParameterSet& p);
		void RecalculateIntegratedHazard(SR::ParameterSet& p);
		void RecalculateAll(SR::ParameterSet& p) {
			RecalculateHazard(p);
			RecalculateLnHazard(p);
			RecalculateIntegratedHazard(p);
		};
		double GetHazard(double t);
		double LnHazard(double t);
		double IntegratedHazard(double t);
	};
	double qsimp(double func(SR::ParameterSet& p, const double), const double a, const double b, SR::ParameterSet& p);
	double trapzd(double func(SR::ParameterSet& p, const double), const double a, const double b, const int n, SR::ParameterSet& p);
};

#endif

