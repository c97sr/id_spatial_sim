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

#include"ebola.h"

extern gsl_rng * glob_rng;

bool evCountSpatialNeighbour(SR::Node* pointer1, SR::Node* pointer2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pointer1->IncrementNoSpatialNeighbour();
	pointer2->IncrementNoSpatialNeighbour();
	return true;
};

bool evCountHouseholdNeighbour(SR::Node* pointer1, SR::Node* pointer2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pointer1->IncrementHouseholdMax();
	pointer2->IncrementHouseholdMax();
	return true;
};

bool evAddToNeighbourList(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pt1->AddToNeighbourList(pt2);
	pt2->AddToNeighbourList(pt1);
	return true;
};

void procSpatialNeighbourSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static SR::UntimedEvent tmpev;
	tmpev.ue = evCountSpatialNeighbour;
	tmpev.ptNode1 = pt1;
	tmpev.ptNode2 = pt2;
	em.AddEvent(tmpev,0);
};

void procSpatialNeighbourAdd(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static SR::UntimedEvent tmpev;
	tmpev.ue = evAddToNeighbourList;
	tmpev.ptNode1 = pt1;
	tmpev.ptNode2 = pt2;
	em.AddEvent(tmpev,0);
};

double kernNeighbourSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* constant = p.GetPointer("Constant_Generate_Spatial_Neighbour");
	static double* decay = p.GetPointer("Decay_Generate_Spatial_Neighbour");
	static double* cutdist = p.GetPointer("Cutoff_Distance_Generate_Spatial_Neighbour");
	static double dist;
	dist = pt1->Distance(pt2) - offset;
	if (dist < *cutdist) return 0;
	return (*constant)*pow((dist+1),-(*decay));
};

double kernIntroSetup(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* prob = p.GetPointer("Prob_Generate_Spatial_Neighbour");
	return (*prob);
};

double kernSpatialInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* constant = p.GetPointer("Relative_Transmit_Spatial");
	static double* decay = p.GetPointer("Decay_Transmit_Spatial");
	static double* cutdist = p.GetPointer("Cutoff_Distance_Transmit_Spatial");
	static double* cuttime = p.GetPointer("Cutoff_Time_Transmit_Spatial");
	static double* cuttimefactor = p.GetPointer("Cutoff_Factor_Transmit_Spatial");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* t = p.GetPointer("dblCurrentTime");
	static double newconstant;
	static double debug;
	double dist;
	newconstant = 1;
	dist = pt1->Distance(pt2)-offset;
	if (*decay == 0) return (*ts)*(*constant);
	if (dist < *cutdist) return 0;
	if (*t < *cuttime) newconstant *= *cuttimefactor;
	debug = (*ts)*newconstant*(*constant)*pow((dist+1),-(*decay));
	return  (*ts)*newconstant*(*constant)*pow((dist+1),-(*decay));
};

double kernSpatialEbola(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* constant = p.GetPointer("Relative_Transmit_Spatial");
	static double* decay = p.GetPointer("Decay_Transmit_Spatial");
	static double* cutdist = p.GetPointer("Cutoff_Distance_Transmit_Spatial");
	static double* cuttime = p.GetPointer("Cutoff_Time_Transmit_Spatial");
	static double* cuttimefactor = p.GetPointer("Cutoff_Factor_Transmit_Spatial");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* t = p.GetPointer("dblCurrentTime");
	static double newconstant;
	static double debug;
	double dist;
	newconstant = 1;
	dist = pt1->Distance(pt2)-offset;
	if (*decay == 0) return (*ts)*(*constant);
	if (dist < *cutdist) return 0;
	if (*t < *cuttime) newconstant *= *cuttimefactor;
	debug = (*ts)*newconstant*(*constant)*pow((dist+1),-(*decay));
	return  (*ts)*newconstant*(*constant)*pow((dist+1),-(*decay));
};

double kernSpatialInfectionOff(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* constant = p.GetPointer("Relative_Transmit_Spatial");
	static double* alpha = p.GetPointer("Decay_Transmit_Spatial");
	static double* rho = p.GetPointer("Offset_Transmit_Spatial");
	static double* cutdist = p.GetPointer("Cutoff_Distance_Transmit_Spatial");
	static double* cuttime = p.GetPointer("Cutoff_Time_Transmit_Spatial");
	static double* cuttimefactor = p.GetPointer("Cutoff_Factor_Transmit_Spatial");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* t = p.GetPointer("dblCurrentTime");
	double dist;
	static double newconstant;
	newconstant = 1;
	if (*alpha == 0) return (*ts)*(*constant);
	dist = pt1->Distance(pt2)-offset;
	if (dist < *cutdist) return 0;
	if (*t < *cuttime) newconstant *= *cuttimefactor;
	return (*ts)*newconstant*(*constant)*(pow((1+dist/(*rho)),-1.0 * (*alpha)));
};

double kernFileCached(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	//  This needs to be the stored convoluted kernel
	static double* constant = p.GetPointer("Relative_Transmit_Spatial");
	static double* h_rrf = p.GetPointer("Hazard_Rash_Relative_Fever");
	static double* ts = p.GetPointer("dblTimeStep");
	double dist;
	static double* cutdist = p.GetPointer("Cutoff_Distance_Transmit_Spatial");
	static double corrected_constant;
	static SR::CachedDblLookup cdl(p.GetTag("strFileMovementKernel"));
	dist = pt1->Distance(pt2);
	dist = dist-offset;
#ifdef _DEBUG
	if (dist < 0) SR::srerror("Negative distance should never be passed to double kernFileCached(SR::Parame...)");
#endif
	if (*cutdist > dist) return 0;
	if (pt1->GetCharacteristic() == 3) corrected_constant = *constant * (*h_rrf);
	else corrected_constant = *constant;
	return  (*ts) * corrected_constant * cdl.GetCheckedValue(dist);
};

double kernFileCachedEarlyRash(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	//  This needs to be the stored convoluted kernel
	static double* constant = p.GetPointer("Relative_Transmit_Spatial");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* rel = p.GetPointer("Hazard_Rash_Relative_Fever");
	double dist;
	static double* cutdist = p.GetPointer("Cutoff_Distance_Transmit_Spatial");
	static SR::CachedDblLookup cdl(p.GetTag("strFileMovementKernel"));
	dist = pt1->Distance(pt2);
	dist = dist-offset;
#ifdef _DEBUG
	if (dist < 0) SR::srerror("Negative distance should never be passed to double kernFileCached(SR::Parame...)");
#endif
	if (*cutdist > dist) return 0;
	return  (*ts) * (*constant) * (*rel) * cdl.GetCheckedValue(dist);
};

double kernHouseholdInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* cc = p.GetPointer("Common_Constant_Transmit");
	static double* hhn = p.GetPointer("Hazard_Home_Relative_Network");
	static double* ts = p.GetPointer("dblTimeStep");
	return 1-exp(-(*cc)*(*hhn)*(*ts));
};

double kernHouseholdInfectionSymp(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* cc = p.GetPointer("Common_Constant_Transmit");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* hrf = p.GetPointer("Hazard_Rash_Relative_Fever");
	static double* hhn = p.GetPointer("Hazard_Home_Relative_Network");
	return 1-exp(-(*cc)*(*hrf)*(*hhn)*(*ts));
};

double kernNeighbourInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* cc = p.GetPointer("Common_Constant_Transmit");
	static double* ts = p.GetPointer("dblTimeStep");
	return 1-exp(-(*cc)*(*ts));
};

double kernNeighbourInfectionSymp(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, double offset) {
	static double* cc = p.GetPointer("Common_Constant_Transmit");
	static double* ts = p.GetPointer("dblTimeStep");
	static double* hrf = p.GetPointer("Hazard_Rash_Relative_Fever");
	return 1-exp(-(*cc)*(*hrf)*(*ts));
};

bool evInfection(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	static int mg=static_cast<int>(p.GetValue("intMaxGeneration"));
	static double* t = p.GetPointer("dblCurrentTime");
	static double* t_mv = p.GetPointer("Start_Time_Movement_Restrictions");
	static double* r_mv = p.GetPointer("Range_Movement_Restrictions");
	static double* e_mv = p.GetPointer("Efficacy_Movement_Restrictions");
	pt2->GetHexagon()->SetMembersAltered();
	if (pt1->GetGeneration() == mg) return false;
	if (pt2->GetCharacteristic() != 0) return false;
	if (pt1->GetQuarantineLevel() > 0) return false;
	if (*t_mv <= *t) {
		if (pt1->Distance(pt2) > 2*(*r_mv)) {
			if (gsl_rng_uniform(glob_rng) < *e_mv) {
//			if (NR::ran2(p.intSeed) < *e_mv) {
				return false;
			}
		}
	}
	pt2->MakeCharacteristicEqualTo(1);
	pt2->SetGeneration(pt1->GetGeneration()+1);
	procSubsequentToInfection(p,pt1,pt2,em);
	return true;
};

bool evVaccinate(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	static double* p_sus_vac_eff = p.GetPointer("Vaccine_Probability_Susceptible");
	static double* lc1_min = p.GetPointer("Latent_Vaccinated_One_Minimum_Time");
	static double* lc1_max = p.GetPointer("Latent_Vaccinated_One_Maximum_Time");
	static double delay;
	static SR::UntimedEvent tmpev;
	pt2->GetHexagon()->SetMembersAltered();
	if (pt2->GetVaccinationClass()==1) return false;
	pt2->SetVaccinationClass(1);
	if (pt2->GetCharacteristic() == 0) {
		if (gsl_rng_uniform(glob_rng) < *p_sus_vac_eff) pt2->MakeCharacteristicEqualTo(5);
// 		if (NR::ran2(p.intSeed) < *p_sus_vac_eff) pt2->MakeCharacteristicEqualTo(5);
		return true;
	}
	if (pt2->GetCharacteristic() == 1) {
		pt2->MakeCharacteristicEqualTo(6);
		delay = 0;
		tmpev.ptNode1 = pt1;
		tmpev.ptNode2 = pt2;
//		delay += static_cast<int>(*lc1_min+static_cast<int>(NR::ran2(p.intSeed)*((*lc1_min-*lc1_max)+1)));
		// delay += static_cast<int>(*lc1_min+static_cast<int>(gsl_rng_uniform(glob_rng)*((*lc1_min-*lc1_max)+1)));
                delay += static_cast<int>(*lc1_min+static_cast<int>(gsl_rng_uniform(glob_rng)*((*lc1_max-*lc1_min)+1)));
		tmpev.ue = evBecomeLatentVaccinatedTwo;
		em.AddEvent(tmpev,static_cast<int>(delay));
		return true;
	}
	if (pt2->GetCharacteristic() == 2) return true;
	return true;
};

bool evBecomeLatentVaccinatedTwo(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	if (pt2->GetCharacteristic() == 6) {
		pt2->MakeCharacteristicEqualTo(7);
		return true;
	}
	return false;
}

bool evBecomeInfectious(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	static double* p_rec_lc1 = p.GetPointer("Recovery_Probability_Latent_Vaccinated_One");
	static double* p_rec_lc2 = p.GetPointer("Recovery_Probability_Latent_Vaccinated_Two");
	if (pt2->GetCharacteristic() == 1) {
		procEnterProdrome(p,pt1,pt2,em);
		return true;
	}
	if (pt2->GetCharacteristic() == 6) {
//		if (NR::ran2(p.intSeed) < *p_rec_lc1) pt2->MakeCharacteristicEqualTo(5);
		if (gsl_rng_uniform(glob_rng) < *p_rec_lc1) pt2->MakeCharacteristicEqualTo(5);
		else procEnterProdrome(p,pt1,pt2,em);
		return true;
	}
	if (pt2->GetCharacteristic() == 7) {
//		if (NR::ran2(p.intSeed) < *p_rec_lc2) pt2->MakeCharacteristicEqualTo(5);
		if (gsl_rng_uniform(glob_rng) < *p_rec_lc2) pt2->MakeCharacteristicEqualTo(5);
		else procEnterProdrome(p,pt1,pt2,em);
		return true;
	}
	if (pt2->GetCharacteristic() == 5) return false;
		SR::srerror("This line should not be reached.\n");
	return false;
};

void procEnterProdrome(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	// Need to edit this so that we get the correct behavior for
	// self isolation onsetting with the start of contact tracing
	static double* pt_mean = p.GetPointer("Fever_Self_Isolation_Average");
	static double* pt_alpha = p.GetPointer("Fever_Self_Isolation_Int_Alpha");
	static double* pt_mean_ar = p.GetPointer("Fever_Self_Isolation_At_Risk_Average");
	static double* pt_alpha_ar = p.GetPointer("Fever_Self_Isolation_Int_At_Risk_Alpha");
	static double* dt = p.GetPointer("dblTimeStep");
	static double* t = p.GetPointer("dblCurrentTime");
	static double* t_si = p.GetPointer("Start_Time_Self_Isolation");
	static double mean,alpha;
	static int max_delay,chosen_delay;
	static SR::UntimedEvent tmpev;
	if (pt2->GetQuarantineLevel()>0) {
		mean = *pt_mean_ar;
		alpha = *pt_alpha_ar;
	} else {
		mean = *pt_mean;
		alpha = *pt_alpha;
	}
	if (*t > *t_si) {
		tmpev.ue = evEnterFeverQuarantine;
		tmpev.ptNode1 = pt1;
		tmpev.ptNode2 = pt2;
		max_delay = em.GetMaxDelay();
		chosen_delay = SR::GammaModelMatrixDelay(mean,static_cast<int>(alpha),p.intSeed,*dt,max_delay-1);
		em.AddEvent(tmpev,chosen_delay);
	}
	pt2->MakeCharacteristicEqualTo(2);
};

bool evEnterEarlyRash(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	static double* pt_mean = p.GetPointer("Rash_Self_Isolation_Average");
	static double* pt_alpha = p.GetPointer("Rash_Self_Isolation_Int_Alpha");
	static double* pt_mean_ar = p.GetPointer("Rash_Self_Isolation_At_Risk_Average");
	static double* pt_alpha_ar = p.GetPointer("Rash_Self_Isolation_Int_At_Risk_Alpha");
	static double* dt = p.GetPointer("dblTimeStep");
	static double* t = p.GetPointer("dblCurrentTime");
	static double* t_si = p.GetPointer("Start_Time_Self_Isolation");
	static double mean,alpha;
	static int max_delay,chosen_delay;
	static SR::UntimedEvent tmpev;
	if (pt2->GetQuarantineLevel()>0) {
		mean = *pt_mean_ar;
		alpha = *pt_alpha_ar;
	} else {
		mean = *pt_mean;
		alpha = *pt_alpha;
	}
	if (*t > *t_si) {
		tmpev.ue = evEnterRashQuarantine;
		tmpev.ptNode1 = pt1;
		tmpev.ptNode2 = pt2;
		max_delay = em.GetMaxDelay();
		chosen_delay = SR::GammaModelMatrixDelay(mean,static_cast<int>(alpha),p.intSeed,*dt,max_delay-1);
		em.AddEvent(tmpev,chosen_delay);
	}
	pt2->MakeCharacteristicEqualTo(3);
	return true;
};

bool evDie(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pt2->MakeCharacteristicEqualTo(4);
	return true;
};

bool evRecover(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pt2->MakeCharacteristicEqualTo(5);
	return true;
};

bool evEndQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pt2->DecrementQuarantineLevel();
	return true;
};

bool evEnterLateRash(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	pt2->MakeCharacteristicEqualTo(8);
	return true;
};

void procInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static int delay;
	static SR::UntimedEvent tmpev;
	delay = 0;
	tmpev.ptNode1 = pt1;
	tmpev.ptNode2 = pt2;
	tmpev.ue = evInfection;
	em.AddEvent(tmpev,delay);
};

void procSubsequentToInfection(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static double* l_mean = p.GetPointer("Latent_Mean");
	static double* l_alpha = p.GetPointer("Latent_Int_Alpha");
	static double* p_mean = p.GetPointer("Prodrome_Mean");
	static double* p_alpha = p.GetPointer("Prodrome_Int_Alpha");
	static double* r_mean = p.GetPointer("Rash_Mean");
	static double* r_alpha = p.GetPointer("Rash_Int_Alpha");
	static double* er_mean = p.GetPointer("Early_Rash_Mean");
	static double* er_alpha = p.GetPointer("Early_Rash_Int_Alpha");
	static double* p_rec = p.GetPointer("Probability_Of_Recovery");
	static double* ts = p.GetPointer("dblTimeStep");
	static int delay;
	static int tmpDelay,tmpDelay2;
	static int maxDelay;
	static SR::UntimedEvent tmpev;
	maxDelay = em.GetMaxDelay();
	delay = 0;
	tmpev.ptNode1 = pt1;
	tmpev.ptNode2 = pt2;

	tmpDelay = SR::GammaModelMatrixDelay(*l_mean,static_cast<int>(*l_alpha),p.intSeed,*ts,maxDelay);
	delay += tmpDelay;
	tmpev.ue = evBecomeInfectious;
	em.AddEvent(tmpev,delay);

	tmpDelay = SR::GammaModelMatrixDelay(*p_mean,static_cast<int>(*p_alpha),p.intSeed,*ts,maxDelay);
	delay += tmpDelay;
	tmpev.ue = evEnterEarlyRash;
	em.AddEvent(tmpev,delay);

	tmpDelay = SR::GammaModelMatrixDelay(*r_mean,static_cast<int>(*r_alpha),p.intSeed,*ts,maxDelay);
	tmpDelay2 = SR::GammaModelMatrixDelay(*er_mean,static_cast<int>(*er_alpha),p.intSeed,*ts,maxDelay);

	delay += tmpDelay;
	// if (NR::ran2(p.intSeed) < *p_rec) tmpev.ue = evDie;
	if (gsl_rng_uniform(glob_rng) < *p_rec) tmpev.ue = evDie;
	else tmpev.ue = evRecover;
	if (delay >= maxDelay) {delay = maxDelay-1;cerr << "Warning: delay too long.\n"; SR::srerror("Delay too long");}
	em.AddEvent(tmpev,delay);

	if (tmpDelay2 < tmpDelay) {
		delay-=tmpDelay;
		delay+=tmpDelay2;
		tmpev.ue = evEnterLateRash;
		em.AddEvent(tmpev,delay);
	}
};

bool evEnterContactQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	static double* p_dq = p.GetPointer("Duration_Quarantine");
	static double* p_dt = p.GetPointer("dblTimeStep");
	static SR::UntimedEvent tmpev2;
	static int delay2;
	delay2=static_cast<int>(*p_dq / *p_dt+0.5);
	pt2->IncrementQuarantineLevel();
	tmpev2.ptNode1=pt1;
	tmpev2.ptNode2=pt2;
	tmpev2.ue = evEndQuarantine;
	em.AddEvent(tmpev2,delay2);
	pt2->GetHexagon()->SetMembersAltered();
	return true;
};

bool evEnterFeverQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	procEnterQuarantine(p,pt1,pt2,em);
	return true;
	/*
	if (pt2->GetCharacteristic()==2) {
		procEnterQuarantine(p,pt1,pt2,em);
		return true;
	} else return false;
	*/
};

bool evEnterRashQuarantine(SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em, SR::ParameterSet& p) {
	procEnterQuarantine(p,pt1,pt2,em);
	return true;
	/*
	if (pt2->GetCharacteristic()==3) {
		procEnterQuarantine(p,pt1,pt2,em);
		return true;
	} else return false;
	*/
};

void procEnterQuarantine(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static double* p_dq = p.GetPointer("Duration_Quarantine");
	static double* p_dt = p.GetPointer("dblTimeStep");
	static double* p_trc_time=p.GetPointer("Start_Time_Contact_Tracing");
	static double* t = p.GetPointer("dblCurrentTime");
	static SR::UntimedEvent tmpev2;
	int delay2;
	delay2 = static_cast<int>(*p_dq / *p_dt+0.5);
	if(*t >= *p_trc_time) procContactTrace(p,pt1,pt2,em); // NMF
	pt2->IncrementQuarantineLevel();
	tmpev2.ptNode1=pt1;
	tmpev2.ptNode2=pt2;
	tmpev2.ue = evEndQuarantine;
	em.AddEvent(tmpev2,delay2);
};

void procContactTrace(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static double* p_house = p.GetPointer("Probability_Trace_Household");
	static double* p_neighbour = p.GetPointer("Probability_Trace_Neighbour");
	static double* p_mean_hh_dv = p.GetPointer("Vaccine_Delay_Household_Average");
	static double* p_alpha_hh_dv = p.GetPointer("Vaccine_Delay_Household_Alpha");
	static double* p_mean_pg_dv = p.GetPointer("Vaccine_Delay_PG_Average");
	static double* p_alpha_pg_dv = p.GetPointer("Vaccine_Delay_PG_Alpha");
	static SR::UntimedEvent tmpev,tmpev2;
	static SR::Node **ptptCurrent,**ptptNeighbourStartHouseFinish,**ptptNeighbourFinish;
	static int delay; // This needs to be drawn froma  random variable
	static int maxDelay;
	maxDelay = em.GetMaxDelay();
	static double* ts = p.GetPointer("dblTimeStep");
	// delay=static_cast<int>(1.0 / *p_dt+0.5);
	tmpev.ptNode1 = pt1;
	tmpev2.ptNode1 = pt1;
	tmpev.ue = evVaccinate;
	tmpev2.ue = evEnterContactQuarantine;
	ptptCurrent = pt2->GetFirstHouseholdMember();
	ptptNeighbourStartHouseFinish = ptptCurrent + pt2->GetHouseholdMax();
	ptptNeighbourFinish = ptptNeighbourStartHouseFinish + pt2->GetNoSpatialNeighbour();
	while (ptptCurrent != ptptNeighbourStartHouseFinish) {
//		if (NR::ran2(p.intSeed) < *p_house) {
		if (gsl_rng_uniform(glob_rng) < *p_house) {
			tmpev.ptNode2 = *ptptCurrent;
			tmpev2.ptNode2 = *ptptCurrent;
			delay = SR::GammaModelMatrixDelay(*p_mean_hh_dv,static_cast<int>(*p_alpha_hh_dv),p.intSeed,*ts,maxDelay);
			em.AddEvent(tmpev,delay);
			em.AddEvent(tmpev2,delay);
		}
		ptptCurrent++;
	}
	while (ptptCurrent != ptptNeighbourFinish) {
//		if (NR::ran2(p.intSeed) < *p_neighbour) {
		if (gsl_rng_uniform(glob_rng) < *p_neighbour) {
			tmpev.ptNode2 = *ptptCurrent;
			tmpev2.ptNode2 = *ptptCurrent;
			delay = SR::GammaModelMatrixDelay(*p_mean_pg_dv,static_cast<int>(*p_alpha_pg_dv),p.intSeed,*ts,maxDelay);
			em.AddEvent(tmpev,delay);
			em.AddEvent(tmpev2,delay);
		}
		ptptCurrent++;
	}
};

void procRegionalVaccination(SR::ParameterSet &p, SR::Node* pt1, SR::Node* pt2, SR::EventMatrix& em) {
	static SR::UntimedEvent tmpev;
	tmpev.ptNode1 = pt1;
	tmpev.ptNode2 = pt2;
	tmpev.ue = evVaccinate;
	em.AddEvent(tmpev,0);
};

string GatherUKPoxInfoA::ConvertEventPointerToString(SR::UNTIMEDEVENT ue) {
	if (ue==evInfection) return "INF";
	if (ue==evDie) return "DTH";
	if (ue==evVaccinate) return "VAC";
	if (ue==evEnterContactQuarantine) return "CTQ";
	return "XXX";
};

GatherUKPoxInfoA::GatherUKPoxInfoA(int nt, double dt, int nr, string s, int ss, int lss) :
vecInc(3,0),
processData(noEvents*nt*nr,0),
vecTimes(nt),
vecAverages(nt),
vecStandardDeviations(nt),
vecCumTotals(nr),
eventStack(ss),
locationStack(ss) {
	dblTimeStep=dt;
	blLogEvents=true;
	maxSpatialExtent=0;
	string s1,s2;
	noRealisations = nr;
	noTimeSteps = nt;
	SetFileBase(s);
	noRealisationsCompleted=0;
	intCurrentTimestep=0;
	intCurrentStackSize=0;
	intCurrentSpatialStackSize=0;
	intLastSpatialStackSize=0;
	s1 = s+"_Events.out";
	s2 = "Run\tDay\tEvent\tIndex\tX\tY\tGeneration\tinfector\tinfect_x\tinfect_y";
	cerr << s1 << "\n";
	if (!SR::OpenNullFile(s1,s2)) SR::srerror("Couldn't initialise event log file in GatherUKPoxInfoA::GatherUKPoxInfoA\n");
};

void GatherUKPoxInfoA::RegisterEventAfterApplication(vector<SR::UntimedEvent>::iterator ptEv, int cb) {
	// Event registration seems to be utrned off??  Should be easily sorted below!
	static vector<double>::iterator ptVal;
	static int index;
	index=-1;
	if (ptEv->ue==evInfection)
		index=0;
#ifndef SR_INFECTIONS_ONLY
	else if (ptEv->ue==evVaccinate)
		index=1;
	else if (ptEv->ue==evDie)
		index=2;
	else if (ptEv->ue==evBecomeInfectious)
		index=4;
	else if (ptEv->ue==evEnterEarlyRash)
		index=5;
	else if (ptEv->ue==evRecover)
		index=6;
	else if (ptEv->ue==evEnterContactQuarantine || ptEv->ue==evEnterFeverQuarantine || ptEv->ue==evEnterRashQuarantine)
		index=7;
#endif // SR_INFECTIONS_ONLY
	else return;
	// Line below increments the event counter
	// Index only has to be consistent with size of process and access data
	// Index =0, 1 or 2 are used for the main events
	// Index = 3 is used for spatial extent elsewhere
	// Index 4 is used for quarantines
	ptVal = AccessData(index,intCurrentTimestep,noRealisationsCompleted);
	(*ptVal)++;
	eventStack[intCurrentStackSize].x=ptEv->ptNode2->GetX();
	eventStack[intCurrentStackSize].y=ptEv->ptNode2->GetY();
	eventStack[intCurrentStackSize].infectx=ptEv->ptNode1->GetX();
	eventStack[intCurrentStackSize].infecty=ptEv->ptNode1->GetY();
	eventStack[intCurrentStackSize].infectorindex=ptEv->ptNode1->GetIndex();
	eventStack[intCurrentStackSize].nodeindex=ptEv->ptNode2->GetIndex();
	eventStack[intCurrentStackSize].realisation=noRealisationsCompleted;
	eventStack[intCurrentStackSize].time=intCurrentTimestep;
	eventStack[intCurrentStackSize].eventindex=index;
	eventStack[intCurrentStackSize].generation=ptEv->ptNode2->GetGeneration();
	intCurrentStackSize++;
	if (intCurrentStackSize==static_cast<int>(eventStack.size())) {
		if (!FlushStack()) SR::srerror("Unable to flush stack in GatherUKPoxInfoA::RegisterEventAfterApplication");
	}
	if (index==0 && AreSpatialMeasuresOn()) {
		AddNodeLocationToSpatialStack(ptEv->ptNode2);
	};
};

bool GatherUKPoxInfoA::WriteAllEventIncidencesToFile() {
	string tmpstr;
	ofstream ofs;
	for (int i=0;i<noEvents;++i) {
		CalcForIndex(i);
		if (i==0) tmpstr = GetFileBase()+"_INF.out";
		if (i==1) tmpstr = GetFileBase()+"_VAC.out";
		if (i==2) tmpstr = GetFileBase()+"_DTH.out";
		if (i==3) tmpstr = GetFileBase()+"_SPX.out";
		if (i==4) tmpstr = GetFileBase()+"_PRO.out";
		if (i==5) tmpstr = GetFileBase()+"_RSH.out";
		if (i==6) tmpstr = GetFileBase()+"_REC.out";
		if (i==7) tmpstr = GetFileBase()+"_CTQ.out";
		ofs.open(tmpstr.c_str());
		if (ofs.fail()) return false;
		ofs << OutputTableForIDL(i);
		ofs.close();
	}
	return true;
}

string GatherUKPoxInfoA::CurrentCumPrev() {
	ostringstream oss;
	for (int i=0;i<static_cast<int>(vecInc.size());++i) oss << vecInc[i] << "\t";
	oss << "\r";
	return oss.str();
};

void GatherUKPoxInfoA::ClearCumPrev() {
	for (int i=0;i<static_cast<int>(vecInc.size());++i) vecInc[i]=0;
};

vector<double>::iterator GatherUKPoxInfoA::AccessData(int ne, int nt, int nr) {
	static int index;
	vector<double>::iterator rtnIt = processData.begin();
	if (ne < 0 || ne >= noEvents)
		SR::srerror("Problem int AccessData A");
	if (nr < 0 || nr >= noRealisations)
		SR::srerror("Problem int AccessData B");
	if (nt < 0 || nt >= noTimeSteps)
		SR::srerror("Problem int AccessData C");
	index = ne + nt*(noEvents) + nr*(noEvents*noTimeSteps);
	if (index >= static_cast<int>(processData.size())) SR::srerror("Index calc not correct.");
	return rtnIt + index;
};

string GatherUKPoxInfoA::OutputRealisationData(int index) {
	ostringstream oss;
	for (int i=0;i<noTimeSteps;++i) {
		for (int j=0;j<noRealisations;++j) oss << *AccessData(index,i,j) << "\t";
		oss << "\n";
	}
	return oss.str();
};

void GatherUKPoxInfoA::UpdateIncidenceFile(int index) {
	ofstream ofs();
};


void GatherUKPoxInfoA::CalcForIndex(int index) {
	static vector<double>::iterator tmp;
	static double cr;
	cr = static_cast<double>(noRealisationsCompleted);
	for (int i=0;i<noTimeSteps;++i) {
		vecAverages[i]=0;
		vecStandardDeviations[i]=0;
		for (int j=0;j<noRealisationsCompleted;++j) {
			tmp = AccessData(index,i,j);
			vecAverages[i]+=*tmp;
			vecStandardDeviations[i]+=(*tmp)*(*tmp);
		}
		vecAverages[i]/=cr;
		vecStandardDeviations[i]=sqrt(vecStandardDeviations[i]/cr-vecAverages[i]*vecAverages[i]);
	}
	for (int i=0;i<noRealisationsCompleted;++i) {
		vecCumTotals[i]=0;
		for (int j=0;j<noTimeSteps;++j) {
			tmp = AccessData(index,j,i);
			vecCumTotals[i]+=*tmp;
		}
	}
	cumAverage=0;
	cumStandardDeviation=0;
	for (int j=0;j<noRealisationsCompleted;++j) {
		cumAverage+=vecCumTotals[j];
		cumStandardDeviation+=vecCumTotals[j]*vecCumTotals[j];
	}
	cumAverage/=cr;
	cumStandardDeviation=sqrt(cumStandardDeviation/cr-cumAverage*cumAverage);
};

string GatherUKPoxInfoA::OutputTable(int index) {
	ostringstream oss;
	oss.precision(5);
	oss << "Compreal\t" << noRealisationsCompleted << "\tnoTimeSteps\t" << noTimeSteps << "\n";
	oss << "Cum\t\t"; for (int i=0;i<noTimeSteps;++i) oss << static_cast<double>(i)*dblTimeStep << "\t"; oss << "\n";
	oss << cumAverage << "\tAve\t"; for (int i=0;i<noTimeSteps;++i) oss << vecAverages[i] << "\t"; oss << "\n";
	oss << cumStandardDeviation << "\tStDev\t"; for (int i=0;i<noTimeSteps;++i) oss << vecStandardDeviations[i] << "\t"; oss << "\n";
	for (int j=0;j<noRealisationsCompleted;++j) {
		oss << vecCumTotals[j] << "\t" << j << "\t";
		for (int i=0;i<noTimeSteps;++i) oss << *(AccessData(index,i,j)) << "\t"; oss << "\n";
	}
	return oss.str();
};

string GatherUKPoxInfoA::OutputTableForIDL(int index) {
	ostringstream oss;
	oss.precision(5);
	oss << noRealisationsCompleted << "\n" << noTimeSteps << "\n";
	for (int i=0;i<noTimeSteps;++i) oss << static_cast<double>(i)*dblTimeStep << "\t"; oss << "\n";
	for (int i=0;i<noTimeSteps;++i) oss << vecAverages[i]/dblTimeStep << "\t"; oss << "\n";
	for (int i=0;i<noTimeSteps;++i) oss << vecStandardDeviations[i]/dblTimeStep << "\t"; oss << "\n";
	for (int j=0;j<noRealisationsCompleted;++j) {
		for (int i=0;i<noTimeSteps;++i) oss << *(AccessData(index,i,j))/dblTimeStep << "\t"; oss << "\n";
	}
	return oss.str();
};

bool GatherUKPoxInfoA::FlushStack() {
	static vector<logEvent>::iterator pt,ptEnd;
	string s = GetFileBase()+"_Events.out";
	if (blLogEvents) {
		ofstream f(s.c_str(),iostream::app);
		// cerr << s << "_END_here\n";
		if (f.fail()) return false;
		f.precision(10);
		pt = eventStack.begin();
		ptEnd = pt + intCurrentStackSize;
		while (pt != ptEnd) {
			f << pt->realisation << "\t" << pt->time*dblTimeStep << "\t" << pt->eventindex << "\t" << pt->nodeindex << "\t"
				<< pt->x << "\t" << pt->y << "\t" << pt->generation << "\t" << pt->infectorindex << "\t" << pt->infectx << "\t" << pt->infecty << "\n";
			pt++;
		}
		f.close();
	}
	intCurrentStackSize=0;
	return true;
};

void GatherUKPoxInfoA::SetFileBase(string s) {
	if (static_cast<int>(s.size())>maxstringlength) SR::srerror("String toio long in GatherUKPoxInfoA::SetFileBase");
	intSizeFileBase = s.size();
	for (int i=0;i<intSizeFileBase;++i) fileBase1[i]=s[i];
};

string GatherUKPoxInfoA::GetFileBase() {
	string rtnstr(fileBase1,intSizeFileBase);
	return rtnstr;
};

bool GatherUKPoxInfoA::AddNodeLocationToSpatialStack(SR::Node *ptN) {
	static vector<coords>::iterator ptC;
	if (intCurrentSpatialStackSize==static_cast<int>(locationStack.size())) return false;
	ptC = locationStack.begin()+intCurrentSpatialStackSize;
	ptC->x = ptN->GetX();
	ptC->y = ptN->GetY();
	intCurrentSpatialStackSize++;
	return true;
};

bool GatherUKPoxInfoA::EventStackFull() {
	if (intCurrentSpatialStackSize==static_cast<int>(locationStack.size())) return true;
	else return false;
};

void GatherUKPoxInfoA::RecalculateSpatialExtent() {
	static double dblTmp;
	static vector<coords>::iterator ptC1,ptC2,ptCStart;
	ptC1 = locationStack.begin()+intLastSpatialStackSize;
	ptC2 = locationStack.begin()+intCurrentSpatialStackSize;
	while (ptC1 != ptC2) {
		ptCStart = locationStack.begin();
		while (ptCStart != ptC1) {
			dblTmp = SR::Distance(ptCStart->x,ptCStart->y,ptC1->x,ptC1->y);
			if (dblTmp>maxSpatialExtent) maxSpatialExtent = dblTmp;
			ptCStart++;
		}
		ptC1++;
	}
	intLastSpatialStackSize = intCurrentSpatialStackSize;
};

void GatherUKPoxInfoA::IncrementTimestep() {
	static vector<double>::iterator ptVal;
	ptVal = AccessData(3,intCurrentTimestep,noRealisationsCompleted);
	*ptVal = maxSpatialExtent;
	RecalculateSpatialExtent();
	intCurrentTimestep++;
};

double DensityTest(double x, double y) {
	static double pi = 3.1415927;
	static double ampx = 0.2;
	static double ampy = 0.0;
	static double base = 1;
	static double xconst = 200.0;
	static double yconst = 1.0;
	static double xoffset = 70;
	static double yoffset = 0;
	double rtnval;
	rtnval = base + ampx*sin(2*pi*(x+xoffset)/xconst)+ ampy*sin(2*pi*(y+yoffset)/yconst);
	if (rtnval < 0) return 0;
	return rtnval;
};

void SetUpExtraParameters(SR::ParameterSet &p) {
	p.ReadParams("Common_Constant_Transmit 1 1"); //This is updated later
	p.ReadParams("intNumberOfBlocks 1 0"); // Calculated later
	p.ReadParams("intBlockSize 1 4096");
	p.ReadParams("intNoCharacteristics 1 9");
	p.ReadParams("intIndexOfMaximalElement 1 0");
	p.ReadParams("intIndexOfSusceptible 1 0");
	p.ReadParams("dblCurrentTime 1 0");
	p.ReadParams("intCurrentRealisation 1 0");
	p.ReadParams("intNumberOfWorkplaces 1 0");
	p.ChangeValue("intNumberOfWorkplaces",static_cast<int>(static_cast<double>(p.GetValue("intNoNodes")) / p.GetValue("dblAverageWorkplaceSize")));
	p.ReadParams("intMaxNoHexagons 1 0");
	//Change these to +1 to fix vecHexagons being too big?
	p.ChangeValue("intMaxNoHexagons", static_cast<int>(p.GetValue("dblXGridSize")/(1.5*p.GetValue("dblHexagonWidth"))+2)*static_cast<int>(p.GetValue("dblYGridSize")/(0.866025*p.GetValue("dblHexagonWidth"))+2));
	p.ReadParams("intNumberOfBlocks 1 0");
	p.ReadParams("dblAveCalcNeighbours 1 0");
	p.ReadParams("dblAveCalcHousehold 1 0");
	p.Lock();
};

void RunExtraParameters(SR::ParameterSet &p) {
	double tmp;
	p.ReadParams("dblStartTime 1 0");
	// These are calculated from the parameters read in from file.
	p.ReadParams("Relative_Transmit_Spatial 1 1");
	p.ReadParams("Common_Constant_Transmit 1 1"); //This is updated later
	p.ReadParams("intMaxDelayTimesteps 1 1");
	tmp = p.GetValue("dblMaxDelayTime")/p.GetValue("dblTimeStep")+1;
	p.ChangeValue("intMaxDelayTimesteps",tmp);
	p.ReadParams("Daily_Expected_Spatial 1 0");
	p.ReadParams("R0_Overall 1 0");
	p.ReadParams("Theta 1 0");
	p.ReadParams("dblPropControlled 1 0");
	p.ReadParams("dblAveTimeThoseControlled 1 0");
	p.ReadParams("dblPropReachedMax 1 0");
	p.ReadParams("dblAveTimeThoseReachedMax 1 0");
	p.ReadParams("dblPropNotControlledNotReachedMax 1 0");
	p.ReadParams("dblAveInfAtEndTimeThoseNotControlledNotReachedMax 1 0");
	p.ReadParams("intAverageInfectionsIfControlled 1 0");
	p.ReadParams("intAverageInfectionsIfNotControlledMaxNotReached 1 0");

	p.Lock();
};

// This function takes the attack rate and tries to calculate the correct attack rates
void CalcBetasWithAttackRate(SR::ParameterSet &p, SR::GridHex& g, SR::KERNEL k) {
	double nn= g.CalculateAverageSpatialNeighbours(); //NMF
	double nh = g.CalculateAverageHousehold();
	cerr << "Ave work contacts (not p*<N_w>) = " << nn << ", ave household contacts (not <N_h>-1) = " << nh << "\n";
	double timestep = p.GetValue("dblTimeStep");
	int maxNoTimesteps = p.GetIntValue("intMaxDelayTimesteps");
	double meanPro = p.GetValue("Prodrome_Mean");
	int alphaPro = p.GetIntValue("Prodrome_Int_Alpha");
	double meanSymp = p.GetValue("Rash_Mean");
	int alphaSymp = p.GetIntValue("Rash_Int_Alpha");
	double meanSympEarly = p.GetValue("Early_Rash_Mean");
	int alphaSympEarly = p.GetIntValue("Early_Rash_Int_Alpha");
	double dblDiscCalcTimestep = 1.0/static_cast<double>(p.GetIntValue("intNoIntStepsPerTimeStep"));
	double h_rf = p.GetValue("Hazard_Rash_Relative_Fever");
	double h_hn = p.GetValue("Hazard_Home_Relative_Network");
	double r0_spatial = p.GetValue("R0_Spatial");
	double r0n = p.GetValue("R0_Network");
	double attack_rate = p.GetValue("AttackRate");
	double ar_accuracy = p.GetValue("dblNetworkAttackRateAccuracy");
	int ar_max_tries = p.GetIntValue("intNetworkAttackRateTries");
	double durProd,durEarlyRash;
	double spatial_coeff;
	int stop_flag=false;
	int notry=0;
	double epsilon=1e-10;
	static double p_symp_home,p_symp_network,p_pro_home,p_pro_network;
	static double lower_beta,current_beta,upper_beta;
	static double current_r0;
	static double lower_ar,current_ar,upper_ar;
	static double current_r0n,lower_r0n,upper_r0n;
	static double theta;

	if (r0n > 0 & attack_rate > 0) SR::srerror("Cannot set R0 network and household attack rate");

	cerr << "Setting vectors of discretised gamma probabilities...";
	vector<double> feverVec = GenerateGammaVector(meanPro,alphaPro,timestep,dblDiscCalcTimestep,maxNoTimesteps);
	vector<double> rashVec = GenerateGammaVector(meanSymp,alphaSymp,timestep,dblDiscCalcTimestep,maxNoTimesteps);
	vector<double> earlyRashVec = GenerateGammaVector(meanSympEarly,alphaSympEarly,timestep,dblDiscCalcTimestep,maxNoTimesteps);
	cerr << "done.\n";

	lower_beta=0;
	upper_beta=100.0;
	lower_ar=0;
	upper_ar=1;

	// Start the different conditions for parameterizing
	// If the attack rate is greater than zero and network R0 zero
	if (attack_rate > epsilon && r0n < epsilon) {
		current_beta=0.5;
		stop_flag=0;
		cerr << "Setting absolute hazards for network transmission using attack rate ...";
		// Should also have the option for setting with R0 here
		while (!stop_flag) {
			p_pro_home = CalcProbContactInfected(feverVec,timestep,current_beta*h_hn);
			p_symp_home = CalcProbContactInfected(rashVec,timestep,current_beta*h_hn*h_rf);
			current_ar = p_pro_home+(1-p_pro_home)*p_symp_home;
			if (attack_rate > current_ar) {lower_ar=current_ar;lower_beta=current_beta;current_beta=(lower_beta+upper_beta)/2.0;}
			else {upper_ar=current_ar;upper_beta=current_beta;current_beta=(lower_beta+upper_beta)/2.0;}
			notry++;
			if (notry > ar_max_tries) SR::srerror("Could not get accurate enough r0 network with allowed number of attempts");
			if (upper_ar - lower_ar < ar_accuracy) {current_beta=(lower_beta+upper_beta)/2.0;stop_flag=true;}
		}
		cerr << "done.\n";

	// Network ro greater than zero but no attack rate
	} else if (attack_rate < epsilon && r0n > epsilon) {
		current_beta=0.5;
		stop_flag=0;
		cerr << "Setting absolute hazards for network transmission using network R0 ...";
		if (nn < epsilon) SR::srerror("Trying to get an r0 in the network without any network!");
		while (!stop_flag) {
			p_pro_network = CalcProbContactInfected(feverVec,timestep,current_beta);
			p_symp_network = CalcProbContactInfected(earlyRashVec,timestep,current_beta*h_rf);
			current_r0n = nn*(p_pro_network+(1-p_pro_network)*p_symp_network);
			if (r0n > current_r0n) {lower_r0n=current_r0n;lower_beta=current_beta;current_beta=(lower_beta+upper_beta)/2.0;}
			else {upper_r0n=current_r0n;upper_beta=current_beta;current_beta=(lower_beta+upper_beta)/2.0;}
			notry++;
			if (notry > ar_max_tries) SR::srerror("Could not get accurate enough r0 network with allowed number of attempts");
			if (upper_r0n - lower_r0n < ar_accuracy*nn) {current_beta=(lower_beta+upper_beta)/2.0;stop_flag=true;}
		}
		cerr << "done.\n";

	// Both r0 and attack rate zero
	} else if (attack_rate < epsilon && r0n < epsilon) {
		current_beta=0;

	// Throw and error if they are both non zero
	} else {
		SR::srerror("Can't have attack rate and network R0 be greater than zero");
	}

	p.ChangeValue("Common_Constant_Transmit",current_beta);
	cerr << "done.\n";

	p_pro_home = CalcProbContactInfected(feverVec,timestep,current_beta*h_hn);
	p_symp_home = CalcProbContactInfected(rashVec,timestep,current_beta*h_hn*h_rf);
	p_pro_network = CalcProbContactInfected(feverVec,timestep,current_beta);
	p_symp_network = CalcProbContactInfected(earlyRashVec,timestep,current_beta*h_rf);

	current_r0 = nh*(p_pro_home+(1-p_pro_home)*p_symp_home)+nn*(p_pro_network+(1-p_pro_network)*p_symp_network);
	current_ar = p_pro_home+(1-p_pro_home)*p_symp_home;

	cerr << "Household attack rate is " <<  current_ar << " with R0 of " << (current_r0+r0_spatial) << ".\n";

	p.ChangeValue("R0_Overall",current_r0+r0_spatial);

	cerr << "Setting absolute hazards for spatial transmission...";

	if (r0_spatial > 0) {
		durProd = SR::MeanOfModelMatrixDelay(meanPro,alphaPro,timestep,dblDiscCalcTimestep,maxNoTimesteps);
		durEarlyRash = SR::MeanOfModelMatrixDelay(meanSympEarly,alphaSympEarly,timestep,dblDiscCalcTimestep,maxNoTimesteps);
		durProd = meanPro;
		durEarlyRash = meanSympEarly;
		spatial_coeff = g.CalcExpectedSpatial(p,k);
			cerr << durEarlyRash << " " << durProd << " " << h_rf << "\n";
		p.ChangeValue("Relative_Transmit_Spatial",1/spatial_coeff);
		p.ChangeValue("Daily_Expected_Spatial",r0_spatial*timestep/(durProd+h_rf*durEarlyRash));
	} else {p.ChangeValue("Relative_Transmit_Spatial",0);}

	theta = (nh*p_pro_home+nn*p_pro_network+durProd*r0_spatial/(durProd+h_rf*durEarlyRash))/(current_r0+r0_spatial);
	p.ChangeValue("Theta",theta);

	cerr << "done.\n";
	cerr << "Beta_network = " << p.GetValue("Common_Constant_Transmit") << " Beta_spatial = " << p.GetValue("Relative_Transmit_Spatial") << "\n";
};

double CalcProbContactInfected(double mean, int alpha, double timestep, double intstep, int maxstep, double hazard) {
	static double prob,currentinfectprob,currentsurviveprob,currentcumsurviveprob,currentoverallprob;
	static int currentstep;
	prob = 1-exp(-hazard*timestep);
	currentinfectprob=0;
	currentsurviveprob=0;
	currentcumsurviveprob=SR::ProbabilityOfOneTimeStep(mean,alpha,timestep,intstep,0);
	currentoverallprob=0;
	currentstep=1;
	while (currentstep <= maxstep) {
		currentinfectprob+=pow(1-prob,currentstep-1)*prob;
		currentsurviveprob=SR::ProbabilityOfOneTimeStep(mean,alpha,timestep,intstep,currentstep);
		currentcumsurviveprob+=currentsurviveprob;
		currentoverallprob+=currentinfectprob*currentsurviveprob;
		currentstep++;
	}
	currentoverallprob+=currentinfectprob*(1-currentcumsurviveprob);
	return currentoverallprob;
};

vector<double> GenerateGammaVector(double mean, int alpha, double timestep, double intstep, int maxstep) {
	vector<double> rtnvec(maxstep+1);
	vector<double>::iterator pt = rtnvec.begin();
	int index = 0;
	while (index <= maxstep) {
		rtnvec[index] = SR::ProbabilityOfOneTimeStep(mean,alpha,timestep,intstep,index);
		index++;
	};
	return rtnvec;
};

double CalcProbContactInfected(vector<double>& gamvec, double timestep, double hazard) {
	static double prob,currentinfectprob,currentsurviveprob,currentcumsurviveprob,currentoverallprob;
	static int currentstep,maxstep;
	maxstep = gamvec.size()-1;
	prob = 1-exp(-hazard*timestep);
	currentinfectprob=0;
	currentsurviveprob=0;
	currentcumsurviveprob=gamvec[0];
	currentoverallprob=0;
	currentstep=1;
	while (currentstep <= maxstep) {
		currentinfectprob+=pow(1-prob,currentstep-1)*prob;
		currentsurviveprob=gamvec[currentstep];
		currentcumsurviveprob+=currentsurviveprob;
		currentoverallprob+=currentinfectprob*currentsurviveprob;
		currentstep++;
	}
	currentoverallprob+=currentinfectprob*(1-currentcumsurviveprob);
	return currentoverallprob;
};

double fnLogLogBiphasic(double d, SR::ParameterSet& p) {
	static double *p1 = p.GetPointer("Commute_Power_One");
	static double *p2 = p.GetPointer("Commute_Power_Two");
	static double *cp = p.GetPointer("Commute_Change_Point");
	if (d < *cp) return pow(d,*p1);
	return pow(*cp,*p1-*p2)*pow(d,*p2);
};

double fnSquaredLogLogBiphasic(double d, SR::ParameterSet& p) {
	static double *p1 = p.GetPointer("Commute_Power_One");
	static double *p2 = p.GetPointer("Commute_Power_Two");
	static double *cp = p.GetPointer("Commute_Change_Point");
	static double debug;
	d = sqrt(d);
	debug = pow(d,*p1);
	if (d < *cp) return pow(d,*p1);
	return pow(*cp,*p1-*p2)*pow(d,*p2);
};

double fnSquaredOffsetPower(double d, SR::ParameterSet& p) {
	static double *p1 = p.GetPointer("Commute_Power_One");
	static double *cp = p.GetPointer("Commute_Change_Point");
	d = sqrt(d);
	return pow(1+d/(*cp),*p1);
};

double fnOffsetPower(double d, SR::ParameterSet& p) {
	static double *p1 = p.GetPointer("Commute_Power_One");
	static double *cp = p.GetPointer("Commute_Change_Point");
	return 1/(1+pow(d/(*cp),*p1));
};

double fnOffsetPower(double d, int np, double *p) {
	static double *p1 = p;		// Commute power
	static double *cp = (p+1);	// Commute change point / offset
	return 1/(1+pow(d/(*cp),*p1));
};

double fnOneOverRadiusOffsetPower(double d, SR::ParameterSet& p) {
	static double *p1 = p.GetPointer("Commute_Power_One");
	static double *cp = p.GetPointer("Commute_Change_Point");
	return 1/d*pow(1+d/(*cp),*p1);
};

double distTestDist(double ave, SR::ParameterSet &p) {
	// switch this to the inverse power law and have the parameter be the variance of the distribution
	static int intAlpha = static_cast<int>(p.GetValue("dblConditionalDistParam"));
	return static_cast<double>(SR::GammaModelDev(ave,intAlpha,p.intSeed));
};

int PoxEventToInt(SR::UNTIMEDEVENT ue) {
	if (ue==evInfection) return 0;
	else if (ue == evVaccinate) return 1;
	else if (ue == evBecomeLatentVaccinatedTwo) return 2;
	else if (ue == evBecomeInfectious) return 3;
	else if (ue == evEnterEarlyRash) return 4;
	else if (ue == evDie) return 5;
	else if (ue == evRecover) return 6;
	else if (ue == evEndQuarantine) return 7;
	else if (ue == evEnterLateRash) return 8;
	else if (ue == evEnterContactQuarantine) return 9;
	else if (ue == evEnterFeverQuarantine) return 10;
	else if (ue == evEnterRashQuarantine) return 11;
	else return 999;
	// else SR::srerror("Event not known here");
};

SR::UNTIMEDEVENT PoxIntToEvent(int e) {
	if (e == 0) return evInfection;
	else if (e == 1) return evVaccinate;
	else if (e == 2) return evBecomeLatentVaccinatedTwo;
	else if (e == 3) return evBecomeInfectious;
	else if (e == 4) return evEnterEarlyRash;
	else if (e == 5) return evDie;
	else if (e == 6) return evRecover;
	else if (e == 7) return evEndQuarantine;
	else if (e == 8) return evEnterLateRash;
	else if (e == 9) return evEnterContactQuarantine;
	else if (e == 10) return evEnterFeverQuarantine;
	else if (e == 11) return evEnterRashQuarantine;
	else SR::srerror("Int not recognized here");
	return evDie;
};

void setGZRegionMembership(SR::GridHex &g) {

	double constA = 3.1415927/180;

	double tmpx[] = {114.487,113.80,113.867,114.027,114.219,114.429};
	double tmpy[] = {22.142,22.151,22.431,22.504,22.547,22.545};
	int novert = 6;

	for (int i=0;i<novert;++i) {
		tmpx[i] = tmpx[i]*constA;
		tmpy[i] = tmpx[i]*constA;
	}

	SR::Node* ptNode = g.FirstNode();
	SR::Node* ptLastNode = ptNode + g.GetNoNodes();

	while (ptNode != ptLastNode) {
		if (SR::pnpoly(4,tmpx,tmpy,ptNode->GetX(),ptNode->GetY())) ptNode->SetKernelIndex(1);
		else ptNode->SetKernelIndex(0);
		ptNode ++;
	}

};
