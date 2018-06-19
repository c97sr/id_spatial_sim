#include"SR_Stats.h"

// Declare global pointer for gsl random numbers
extern gsl_rng * glob_rng;

double SR::LogNormalPdf(double *x, double *mu, double *sigma) {
	static double rtnval;
	if (*x < 1e-100) return 0;
	rtnval = 1/(2.506628274631*(*sigma)*(*x))*exp(-(log(*x)-(*mu))*(log(*x)-(*mu))/(2*(*sigma)*(*sigma)));
	return rtnval;
	// return 1.0/(*x+1.0);
};

double SR::LnGamma(double *x, double *lambda, double *alpha) {
	// Taken from p50 of Rice, Mathematical Statistics and Data Analysis, 2nd Ed
	static double rtnval;
//	rtnval = (*alpha)*log(*lambda)-NR::gammln(*alpha)+(*alpha - 1)*log(*x) - (*lambda)*(*x);
	rtnval = (*alpha)*log(*lambda)-gsl_sf_lngamma(*alpha)+(*alpha - 1)*log(*x) - (*lambda)*(*x);
	return rtnval;
};

double SR::GammaExcel(double x, double alpha, double beta) {
	// parameterization from MS excel (see "More Help" for GAMMADIST function in Excel)
	static double rtnval;
//	rtnval = 1/pow(beta,alpha)/exp(NR::gammln(alpha))*pow(x,(alpha-1))*exp(-x/beta);
	rtnval = 1/pow(beta,alpha)/exp(gsl_sf_lngamma(alpha))*pow(x,(alpha-1))*exp(-x/beta);
	return rtnval;
};

double SR::GammaModel(double x, double mean, double alpha) {
	// My preferred parameterization.  Overall mean and number of exponentials required (if integer alpha)
	// remember var = mean*mean/alpha
	static double beta;
	beta = mean/alpha;
	return SR::GammaExcel(x,alpha,beta);
};

double SR::GammaModelDev(double mean, int alpha, int &sd) {
	// Random deviate from gamma distribution GammaModel (with integer alphas)
	// remember var = mean*mean/alpha
	static double rtnval;
	static double beta;
	// static int counter = 0;
	// static int sumdev = 0;
	// static int sumsqu = 0;

	beta = mean/static_cast<double>(alpha);

	// rtnval=NR::gamdev(alpha,sd);
	// rtnval*=mean/alpha;

	// XX This one needs a short test case to make sure its working
	rtnval=gsl_ran_gamma(glob_rng,static_cast<double>(alpha),beta);

	// These were the lines used to check the swap out of NRcpp
	// counter = counter + 1;
	// sumdev = sumdev + rtnval;
	// sumsqu = sumsqu + rtnval * rtnval;
	// cerr 	<< counter << "\t" << static_cast<double>(sumdev)/static_cast<double>(counter) << " "
	// 		<< static_cast<double>(sumsqu)/static_cast<double>(counter) << "\n";

	return rtnval;

};

int SR::GammaModelMatrixDelay(double mean, int alpha, int &sd, double dt, int max){
	static int wait;
	static double mean_dash;
	mean_dash = mean/dt+0.5;
	wait = static_cast<int>(SR::GammaModelDev(mean_dash,alpha,sd));
	// wait = static_cast<int>(SR::GammaModelDev(mean,alpha,sd));
	// wait = static_cast<int>(mean/dt);
	// cerr << mean << "\t" << rtnval*dt << "\n";  // this seems OK
	// wait = NR::poidev(mean_dash,sd);
	// if (wait > max) wait = max;
	if (wait > max-1) {
		SR::srerror("Max length exceeded in discrete distribution.");
	}
	return wait;
};

int SR::GammaModelMatrixDelayFixed(double mean, int alpha, int &sd, double dt, int max){
	static int wait;
	wait = static_cast<int>(mean/dt);
	return wait;
};

double SR::MeanOfModelMatrixDelay(double mean, int alpha, double dt_disc, double dt_int, int maxSteps) {
	// routine to give the equivalent mean
	cerr << "Calculating actual means for discrete gamma distribution...";
	static double current_p,running_p,rtnmean,next_dt;
	static int intSteps=0;
	running_p=0;
	rtnmean=0;
	intSteps = 0;
	next_dt = 0;
	while (intSteps <= maxSteps) {
		next_dt += dt_disc;
		current_p = SR::ProbabilityOfOneTimeStep(mean,alpha,dt_disc,dt_int,intSteps);
		running_p += current_p;
		rtnmean += current_p*(next_dt-dt_disc);
		intSteps++;
	}
	rtnmean+=(1-running_p)*(next_dt-dt_disc);
	cerr << "done.\n";
	return rtnmean;
};

double SR::ProbabilityOfOneTimeStep(double mean, int alpha, double dt_disc, double dt_int, int timeStep) {
	static double rtnval,pdf;
	static double current_time;
	static double end_time;
	static double beta;
	static double t;
	static double mean_dash;
	rtnval=0;
	mean_dash = mean + dt_disc/2;
	beta = mean_dash / alpha;
	current_time = static_cast<double>(timeStep)*dt_disc;
	end_time = current_time + dt_disc;
	while (current_time <= end_time) {
		t = current_time + dt_int/2;
//		pdf = 1/pow(beta,alpha)/exp(NR::gammln(alpha))*pow(t,(alpha-1))*exp(-t/beta);
		pdf = 1/pow(beta,alpha)/exp(gsl_sf_lngamma(alpha))*pow(t,(alpha-1))*exp(-t/beta);
		rtnval += pdf*dt_int;
		current_time += dt_int;
	}
	return rtnval;
};

double SR::ignbin(double pp,int nint, int& sd){
	// Loads of comments here
	static double psave = -1.0E37;
	static double nsave = -214748365;
	static double i,k,mp,T1,ix,ix1,m;
	static double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
		x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;
	double n;

	n=(double) nint;
	if(pp != psave) goto S10;
	if(n != nsave) goto S20;
	if(xnp < 30.0) goto S150;
	goto S30;
S10:
	// *****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
	// JJV added checks to ensure 0.0 <= PP <= 1.0
	if(pp < 0.0) SR::srerror("PP < 0.0 in IGNBIN");
	if(pp > 1.0) SR::srerror("PP > 1.0 in IGNBIN");
	psave = pp;
	p = ((psave) <= (1.0-psave) ? (psave) : (1.0-psave));
	q = 1.0-p;
S20:
	// JJV added check to ensure N >= 0
	if(n < 0)
		SR::srerror("N < 0 in IGNBIN");
	xnp = n*p;
	nsave = n;
	if(xnp < 30.0) goto S140;
	ffm = xnp+p;
	m = floor(ffm);
	fm = m;
	xnpq = xnp*q;
	p1 = (int) (2.195*sqrt(xnpq)-4.6*q)+0.5;
	xm = fm+0.5;
	xl = xm-p1;
	xr = xm+p1;
	c = 0.134+20.5/(15.3+fm);
	al = (ffm-xl)/(ffm-xl*p);
	xll = al*(1.0+0.5*al);
	al = (xr-ffm)/(xr*q);
	xlr = al*(1.0+0.5*al);
	p2 = p1*(1.0+c+c);
	p3 = p2+c/xll;
	p4 = p3+c/xlr;
S30:
	// *****GENERATE VARIATE
	u = NR::ran2(sd)*p4;
	v = NR::ran2(sd);
	//TRIANGULAR REGION
	if(u > p1) goto S40;
	ix = floor(xm-p1*v+u);
	goto S170;
S40:

	//     PARALLELOGRAM REGION

	if(u > p2) goto S50;
	x = xl+(u-p1)/c;
	v = v*c+1.0-((xm-x)>= 0?(xm-x):-(xm-x))/p1;
	if(v > 1.0 || v <= 0.0) goto S30;
	ix = floor(x);
	goto S70;
S50:

	//     LEFT TAIL

	if(u > p3) goto S60;
	ix = floor(xl+log(v)/xll);
	if(ix < 0) goto S30;
	v *= ((u-p2)*xll);
	goto S70;
S60:

	//     RIGHT TAIL

	ix = floor(xr-log(v)/xlr);
	if(ix > n) goto S30;
	v *= ((u-p3)*xlr);
S70:

	//*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST

	k =fabs(ix-m);
	if(k > 20 && k < xnpq/2-1) goto S130;

	//     EXPLICIT EVALUATION

	f = 1.0;
	r = p/q;
	g = (n+1)*r;
	T1 = m-ix;
	if(T1 < 0) goto S80;
	else if(T1 == 0) goto S120;
	else  goto S100;
S80:
	mp = m+1;
	for(i=mp; i<=ix; i++) f *= (g/i-r);
	goto S120;
S100:
	ix1 = ix+1;
	for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
	if(v <= f) goto S170;
	goto S30;
S130:

	//     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))

	amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
	ynorm = -(k*k/(2.0*xnpq));
	alv = log(v);
	if(alv < ynorm-amaxp) goto S170;
	if(alv > ynorm+amaxp) goto S30;

	//     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
	//     THE FINAL ACCEPTANCE/REJECTION TEST

	x1 = ix+1.0;
	f1 = fm+1.0;
	z = n+1.0-fm;
	w = n-ix+1.0;
	z2 = z*z;
	x2 = x1*x1;
	f2 = f1*f1;
	w2 = w*w;
	if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
		(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
		(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
		(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
		-140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
	goto S30;
S140:

	//     INVERSE CDF LOGIC FOR MEAN LESS THAN 30

	qn = pow(q,(double)n);
	r = p/q;
	g = r*(n+1);
S150:
	ix = 0;
	f = qn;
	u = NR::ran2(sd);
S160:
	if(u < f) goto S170;
	if(ix > 110) goto S150;
	u -= f;
	ix += 1;
	f *= (g/ix-r);
	goto S160;
S170:
	if(psave > 0.5) ix = n-ix;
	if(floor(ix)!=ix) fprintf(stderr,"*** Non int ignbin ***\n");
	return ix;
};

double SR::rngtest() {

	// const gsl_rng_type * T;

    // Setup random number generator before ever used
	// gsl_rng_env_setup();
	// T = gsl_rng_default;
	// glob_rng = gsl_rng_alloc (T);

	int i, n = 5;

	for (i = 0; i < n; i++) {
		double u = gsl_rng_uniform (glob_rng);
		cout << u << endl;
	}

	return 2.3;

};

