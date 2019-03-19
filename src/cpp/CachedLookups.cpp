/*  Copyright 2018 Steven Riley.

    This file is part of id_spatial_sim.

    id_spatial_sim is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    id_spatial_sim is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with id_spatial_sim.  If not, see <https://www.gnu.org/licenses/>.  */

#include"CachedLookups.h"

SR::CachedIntLookup::~CachedIntLookup() {
	delete [] values;
}

SR::CachedIntLookup::CachedIntLookup(double f(int,SR::ParameterSet&),int min_in,int max_in, SR::ParameterSet& p) {
	values = new double[max_in-min_in+1];
	if (max_in < min_in) SR::srerror("Max can't be less than min in SR::CachedIntLookup");
	max = max_in;
	min = min_in;
	ptStart = values;
	double* ptValue = ptStart;
	int current = min;
	while (current != max+1) {
		*ptValue = f(current,p);
		current++; ptValue++;
	}
	if (ptValue != values+intNoValues) SR::srerror("Problem in SR::CachedIntLookup");
};

double SR::CachedIntLookup::GetCheckedValue(int i) {
	if ((i < min) || (i > max))
		SR::srerror("Out of range in double SR::CachedIntLookup::GetCheckedValue(int i)");
	return *(ptStart+i);
};

double SR::CachedIntLookup::GetFastValue(int i) {
	return *(ptStart+i);
};

SR::CachedDblLookup::CachedDblLookup(double f(double,SR::ParameterSet&),double start_in, double incs_per_order, double no_orders, SR::ParameterSet& p) {
	values = new double[max_values];
	BeenMadeHistogram=false;
	values_used = static_cast<int>(incs_per_order*no_orders+1);
	start=start_in; logstart = log10(start);
	if (values_used >= max_values || values_used < 10) SR::srerror("Too many or too few increments used in CachedDblLookup::CachedDblLookup(double f(double,SR...");
	dlnx = 1/incs_per_order;
	double log_value,linear_value;
	double* ptValue = values;
	*ptValue = f(0,p);ptValue++;
	log_value = log10(start);
	while (ptValue != values+values_used) {
		linear_value = pow(10.0,log_value);
		*ptValue = f(linear_value,p);
		log_value+=dlnx;
		ptValue++;
	}
	max = pow(10.0,log_value-dlnx);
};

SR::CachedDblLookup::CachedDblLookup(double f(double,int, double*),double start_in, double incs_per_order, double no_orders,int np, double *vecp) {
	values = new double[max_values];
	BeenMadeHistogram=false;
	values_used = static_cast<int>(incs_per_order*no_orders+1);
	start=start_in; logstart = log10(start);
	if (values_used >= max_values || values_used < 10) SR::srerror("Too many or too few increments used in CachedDblLookup::CachedDblLookup(double f(double,SR...");
	dlnx = 1/incs_per_order;
	double log_value,linear_value;
	double* ptValue = values;
	*ptValue = f(0,np,vecp);ptValue++;
	log_value = log10(start);
	while (ptValue != values+values_used) {
		linear_value = pow(10.0,log_value);
		*ptValue = f(linear_value,np,vecp);
		log_value+=dlnx;
		ptValue++;
	}
	max = pow(10.0,log_value-dlnx);
};

SR::CachedDblLookup::~CachedDblLookup() {
	delete [] values;
}

SR::CachedDblLookup::CachedDblLookup(string filename) {
	values = new double[max_values];
	string junk;
	double incs_per_order;
	double *ptCurrent = values;
	ifstream ifs(filename.c_str());
	if (ifs.fail()) SR::srerror("Problem opening file in SR::CachedDblLookup::CachedDblLookup(string filename,SR::ParameterSet& p)");
	ifs >> junk >> start;
	logstart=log10(start);
	ifs >> junk >> incs_per_order;
	dlnx = 1/incs_per_order;
	max=logstart;
	ifs >> *ptCurrent;
	values_used=1;ptCurrent++;
	while (ifs >> *ptCurrent) {
		values_used++;
		ptCurrent++;
		max+=dlnx;
		if (ptCurrent == values + max_values) SR::srerror("Not enough values in file in SR::CachedDblLookup::CachedDblLookup(string filename,SR::ParameterSet& p)");
	}
	max = pow(10.0,max);
	ifs.close();
};

double SR::CachedDblLookup::GetCheckedValue(double d) {
	static int index;
	static double lb,ub,lx,logd,rtnval;
	static double epsilon = 1e-10;
	if ( d < epsilon )
		SR::srerror("Problem in SR::CachedDblLookup::GetCheckedValue(double d)");
	if (d > max) {
		cerr << "Warning: cached kernel distance too long.";
		d=max;
	}
	logd = log10(d);
	if (logd < logstart) {
		lb = values[0];
		ub = values[1];
		return lb + (ub-lb)*d;
	}
	index = static_cast<int>((logd-logstart)/dlnx)+1;
	lb = values[index];
	ub = values[index+1];
	lx = logstart+(index-1)*dlnx;
	rtnval = lb + (ub-lb)*(logd-lx)/dlnx;
	return rtnval;
};

double SR::CachedDblLookup::GetHistValue(double d) {
	static int index;
	static double logd;
	static ostringstream oss;
	if (!BeenMadeHistogram)
		SR::srerror("Needs to have been made a histogram in SR::CachedDblLookup::GetHistValue(double d)");
	if ( d >= max || d < start ) {
		oss.str("");
		oss << d;
		SR::srerror("Problem in SR::CachedDblLookup::GetHistValue(double d): "+oss.str());
	}
	logd = log10(d);
	index = static_cast<int>((logd-logstart)/dlnx);
	return values[index];
};

void SR::CachedDblLookup::MakeHistogram() {
	double lnxvalue = log10(start);
	double lineargap;
	for (int i=0;i<values_used-1;++i) {
		lineargap = pow(10.0,lnxvalue+dlnx)-pow(10.0,lnxvalue);
		values[i]=(values[i]+values[i+1])/2*lineargap;
	}
	values[values_used]=0;
	BeenMadeHistogram=true;
};

void SR::CachedDblLookup::Normalise() {
	double total_prob=0;
	for (int i=0;i<values_used;++i) total_prob+=values[i];
	for (int i=0;i<values_used;++i) values[i]/=total_prob;
};

string SR::CachedDblLookup::PrintDistributionOfCommutes() {
	ostringstream oss;
	for (int i=0;i<values_used;++i) {
		oss << pow(10.0,log10(start)+1.0*i*dlnx) << "\t" << values[i] << "\n";
	}
	return oss.str();
};
