#include"SR_Hazard.h"
#include <cmath>
#include "nr.h"

using namespace std;

double SR::qsimp(double func(SR::ParameterSet& p, const double), const double a, const double b, SR::ParameterSet& p)
{
	const int JMAX=20;
	const double EPS=1.0e-10;
	int j;
	double s,st,ost=0.0,os=0.0;

	for (j=0;j<JMAX;j++) {
		st=trapzd(func,a,b,j+1,p);
		s=(4.0*st-ost)/3.0;
		if (j > 5)
			if (fabs(s-os) < EPS*fabs(os) ||
				(s == 0.0 && os == 0.0)) return s;
		os=s;
		ost=st;
	}
	SR::srerror("Too many steps in routine qsimp");
	return 0.0;
}

double SR::trapzd(double func(SR::ParameterSet& p, const double), const double a, const double b, const int n, SR::ParameterSet& p)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(func(p,a)+func(p,b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(p,x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

SR::Hazard::Hazard(double (*h_in)(SR::ParameterSet&, double), int resolution, double min_in, double max_in, SR::ParameterSet& p) {
	dblMin = min_in;
	dblMax = max_in;
	if (dblMax <= dblMin) SR::srerror("Max must be greater than min in Hazard");
	dblDT = (dblMax-dblMin)/static_cast<double>(resolution);
	intNoPoints = resolution+2;
	ptFirstHazard = new double[intNoPoints];
	ptFirstLnHazard = new double[intNoPoints];
	ptFirstIntegratedHazard = new double[intNoPoints];
	h=h_in;
	RecalculateHazard(p);
	RecalculateLnHazard(p);
	RecalculateIntegratedHazard(p);
};

SR::Hazard::~Hazard(){
	delete [] ptFirstHazard;
	delete [] ptFirstLnHazard;
	delete [] ptFirstIntegratedHazard;
};		

void SR::Hazard::RecalculateHazard(SR::ParameterSet& p) {
	static double* ptDbl;
	ptDbl = ptFirstHazard;
	for (double i=0;i<intNoPoints;i++) {
		*ptDbl = h(p,dblMin+i*dblDT);
		ptDbl++;
	}
};

void SR::Hazard::RecalculateLnHazard(SR::ParameterSet& p) {
	static double* ptDbl;
	ptDbl = ptFirstLnHazard;
	for (double i=0;i<intNoPoints;i++) {
		*ptDbl = log(h(p,dblMin+i*dblDT));
		ptDbl++;
	}
};

void SR::Hazard::RecalculateIntegratedHazard(SR::ParameterSet& p) {
	static double* ptDbl;
	static double cumulative;
	cumulative = 0;
	ptDbl = ptFirstIntegratedHazard;
	*ptDbl = cumulative;
	ptDbl++;
	for (double i=1;i<intNoPoints;i++) {
		cumulative += qsimp(h,dblMin+(i-1)*dblDT,dblMin+i*dblDT,p);
		*ptDbl = cumulative;
		ptDbl++;
	}
};

double SR::Hazard::GetHazard(double t) {
	return LinearInterpolate(ptFirstHazard,t);
};

double SR::Hazard::LnHazard(double t) {
	return LinearInterpolate(ptFirstLnHazard,t);
};

double SR::Hazard::IntegratedHazard(double t) {
	return LinearInterpolate(ptFirstIntegratedHazard,t);
};

double SR::Hazard::LinearInterpolate(double *ptM, double t) {
	static double remainder,difference,rtnval,tmpDbl;
	static int index;
#ifdef _DEBUG 
	if (t<0 || t>dblMax) SR::srerror("Bounds problem in SR::Hazard::LnHazard(double t)");
#endif
	if (t<dblMin) t = dblMin;
	index = static_cast<int>((t-dblMin)/dblDT);
	remainder = t-static_cast<double>(index)*dblDT;
	tmpDbl = ptM[index];
	difference = ptM[index+1]-tmpDbl;
	rtnval = tmpDbl+remainder/dblDT*difference;
	return rtnval;
};

