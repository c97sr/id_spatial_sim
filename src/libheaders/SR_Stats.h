#ifndef SR_INC_STATS
#define SR_INC_STATS

#include<iostream>
#include"nr.h"
#include"SR_Utility.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_randist.h>

/*
[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {

	double LnGamma(double *x, double *lambda, double *alpha);
	double LogNormalPdf(double *x, double *mu, double *sigma);
	double ignbin(double pp,int nint, int& sd);
	double GammaExcel(double x, double alpha, double beta);
	double GammaModel(double x, double mean, double alpha);
	double GammaModelDev(double mean, int alpha, int &sd);
	int GammaModelMatrixDelay(double mean, int alpha, int &sd, double dt, int max);
	int GammaModelMatrixDelayFixed(double mean, int alpha, int &sd, double dt, int max);
	double MeanOfModelMatrixDelay(double mean, int alpha, double dt_disc, double dt_int, int max);
	double ProbabilityOfOneTimeStep(double mean, int alpha, double dt_disc, double dt_int, int timeStep);
	double rngtest();

	class GslWrapper {
	private:
		gsl_rng * glob_rng;
		const gsl_rng_type * T;
	public:
		GslWrapper(int seed);
		double GetUniform();
	};

}

#endif
