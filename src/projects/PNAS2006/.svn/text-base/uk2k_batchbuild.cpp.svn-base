#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<fstream>
#include<cmath>
#include<algorithm>

using namespace std;

namespace SR {
	class TemporalData {
	private:
		vector<double> rawData;
		int noRealisations,noTimeSteps;
	public:
		TemporalData(){cerr<<"No defaul constructor for SR::TemporalData";exit(1);};
		TemporalData(int maxNoReal, int maxNoSteps):rawData((maxNoReal+2)*(maxNoSteps+1),-999) {};
		void LoadFromFile(string filename);
		double GetLastAverage();
		pair<double, double> GetLastPercentiles(double lp, double up);
	};
}

int main(int argc, char* argv[]) {

	if (argc!=2) {cerr << "One argmuent required"; exit(1);}
	string arg=argv[1];

	string CMD = "bsub ";
	string srjHOME = "\\\\dide\\usr\\simulation\\Bioterrorism";
	string srjNETS = srjHOME + "\\networks";
	string srjBATCH = srjHOME + "\\batchfiles";
	string srjBIN = srjHOME + "\\bin";
	string srjRES = srjHOME + "\\results";
	string srjPAR = srjHOME + "\\parameters";
	string srjDATA = srjHOME + "\\data";

	SR::TemporalData srtp(1000,1000);
	// srtp.LoadFromFile(srjRES+"\\"+"luton2_INF.txt");

	vector<string> vpPD,vpPDn,vpSD,vpWS,vpGS,vpR0,vpMG,vpMT,vpRP,vpPS,vpTP,vpHP,vpHS;
	ostringstream ossf,osso;
	ofstream ofs1,ofs2;
	string batchname;

	if (arg=="SETUP") {
	cerr << "Making setup script\n";
	vpPD.resize(0);
	vpPDn.resize(0);
	vpSD.resize(0);
	vpWS.resize(0);
	vpGS.resize(0);
	vpR0.resize(0);
	vpMG.resize(0);
	vpMT.resize(0);
	vpRP.resize(0);
	vpPS.resize(0);
	vpTP.resize(0);
	vpHP.resize(0);
	vpHS.resize(0);

	vpPD.push_back(srjDATA+"\\EWPopDense.txt "+srjDATA +"\\EWWorkDense.txt ");
	vpPD.push_back(srjDATA+"\\FlatDense.txt "+srjDATA +"\\FlatDense.txt ");
	vpPDn.push_back("RD");
	vpPDn.push_back("FD");
	vpSD.push_back("3984029522");
	vpSD.push_back("4958948958");
	vpSD.push_back("2309430938");
	vpSD.push_back("5834729388");
	vpSD.push_back("1928499182");
	vpSD.push_back("0594900321");
	vpWS.push_back("5");
	vpWS.push_back("10");
	vpWS.push_back("15");
	vpWS.push_back("20");
	vpGS.push_back("6.392809287");
	vpGS.push_back("14.29475613");
	vpGS.push_back("20.21583799");
	vpGS.push_back("45.20398798");
	vpGS.push_back("63.92809287");
	vpR0.push_back("4");
	vpMG.push_back("4");
	vpMT.push_back("100");
	vpRP.push_back("5");
	vpPS.push_back("1");
	vpTP.push_back("0.2");
	vpHP.push_back("0.25");
	vpHS.push_back("0.125");

	batchname = srjBATCH+"\\"+"renameme_setup.bat";
	ofs1.open(batchname.c_str());
	if (ofs1.fail()) {cerr << "Error opening output file " << batchname << "\n";exit(1);};
	for (int h=0;h<vpPD.size();++h) {
		for (int i=0;i<vpGS.size();++i) {
			for (int j=0;j<vpWS.size();++j) {
				for (int k=0;k<vpSD.size();++k) {
					for (int l=0;l<vpR0.size();++l) {
						for (int m=0;m<vpMG.size();++m) {
							for (int n=0;n<vpMT.size();++n) {
								for (int o=0;o<vpRP.size();++o) {
									for (int p=0;p<vpPS.size();++p) {
										for (int q=0;q<vpTP.size();++q) {
											for (int r=0;r<vpHP.size();++r) {
												for (int s=0;s<vpHS.size();++s) {
													ossf.str("");
													ossf << vpPDn[h] << "_GS_" << vpGS[i] << "_WS_" << vpWS[j] << "_SD_" << vpSD[k];
													ofs1  << CMD
														 << srjBIN << "\\uk2k_build.exe " 
														 << srjPAR << "\\paramBasic.txt " 
														 << vpPD[h] << " "
														 << srjNETS << "\\" << ossf.str() << " "
														 << "Seed " << vpSD[k] << " "
														 << "dblAverageWorkplaceSize " << vpWS[j] << " "
														 << "dblXGridSize" << " " << vpGS[i] << " "
														 << "dblYGridSize" << " " << vpGS[i] << " "
														 << "dblXGridMin 500 dblYGridMin 194 "
														 << "R0_Total " << vpR0[l] << " "
														 << "intMaxGeneration " << vpMG[m] << " "
														 << "EndTime " << vpMT[n] << " "
														 << "RealisationsPerParameterSet " << vpRP[o] << " "
														 << "R0_Proportion_Spatial " << vpPS[p] << " "
														 << "Relative_Transmit_Prodromal " << vpPS[q] << " "
 														 << "Relative_Transmit_Neighbour_vs_House_Prodromal " << vpHP[r] << " "
 														 << "Relative_Transmit_Neighbour_vs_House_Symptomatic " << vpHS[s] << " "
														 << "\n";
												}
											}
										}
									}
								}
							}
						}					
					}
				}
			}
		}
	}
	ofs1 << "pause\n";
	ofs1.close();
	}

	if (arg=="FIG2") {
	cerr << "Making batch script for fig 2 results\n";
	vpPD.resize(0);
	vpPDn.resize(0);
	vpSD.resize(0);
	vpWS.resize(0);
	vpGS.resize(0);
	vpR0.resize(0);
	vpMG.resize(0);
	vpMT.resize(0);
	vpRP.resize(0);
	vpPS.resize(0);
	vpTP.resize(0);
	vpHP.resize(0);
	vpHS.resize(0);

	vpPD.push_back(srjDATA+"\\EWPopDense.txt "+srjDATA +"\\EWWorkDense.txt ");
	vpPD.push_back(srjDATA+"\\FlatDense.txt "+srjDATA +"\\FlatDense.txt ");
	vpPDn.push_back("RD");
	vpPDn.push_back("FD");
	vpSD.push_back("3984029522");
	// vpSD.push_back("4958948958");
	// vpSD.push_back("2309430938");
	// vpSD.push_back("5834729388");
	// vpSD.push_back("1928499182");
	// vpSD.push_back("0594900321");
	vpWS.push_back("5");
	// vpWS.push_back("10");
	// vpWS.push_back("15");
	vpWS.push_back("20");
	vpGS.push_back("6.392809287");
	vpGS.push_back("14.29475613");
	vpGS.push_back("20.21583799");
	vpGS.push_back("45.20398798");
	vpGS.push_back("63.92809287");
	vpR0.push_back("4");
	vpR0.push_back("6");
	vpR0.push_back("8");
	vpMG.push_back("4");
	vpMT.push_back("100");
	vpRP.push_back("1000");
	vpPS.push_back("1");
	vpPS.push_back("0.1");
	vpTP.push_back("0.2");
	vpHP.push_back("0.25");
	vpHS.push_back("0.125");

	batchname = srjBATCH+"\\"+"renameme_fig2.bat";
	ofs2.open(batchname.c_str());
	if (ofs2.fail()) {cerr << "Error opening output file " << batchname << "\n";exit(1);};
	for (int h=0;h<vpPD.size();++h) {
		for (int i=0;i<vpGS.size();++i) {
			for (int j=0;j<vpWS.size();++j) {
				for (int k=0;k<vpSD.size();++k) {
					for (int l=0;l<vpR0.size();++l) {
						for (int m=0;m<vpMG.size();++m) {
							for (int n=0;n<vpMT.size();++n) {
								for (int o=0;o<vpRP.size();++o) {
									for (int p=0;p<vpPS.size();++p) {
										for (int q=0;q<vpTP.size();++q) {
											for (int r=0;r<vpHP.size();++r) {
												for (int s=0;s<vpHS.size();++s) {
													ossf.str("");
													ossf << vpPDn[h] << "_GS_" << vpGS[i] << "_WS_" << vpWS[j] << "_SD_" << vpSD[k] << ".hex";
													osso.str("");
													osso << ossf.str() << "_R0_" << vpR0[l] << "_MG_" << vpMG[m] << "_MT_" << vpMT[n] << "_RP_" << vpRP[o] << "_PS_" << vpPS[p] << "_TP_" << vpTP[q] << "_HP_" << vpHP[r] << "_HS_" << vpHS[s];
													ofs2  << CMD
														 << srjBIN << "\\uk2k_run.exe " 
														 << srjPAR << "\\paramBasic.txt " 
														 << srjNETS << "\\" << ossf.str() << " "
														 << srjRES << "\\" << osso.str() << " "
														 << "Seed " << vpSD[k] << " "
														 << "dblAverageWorkplaceSize " << vpWS[j] << " "
														 << "dblXGridSize" << " " << vpGS[k] << " "
														 << "dblYGridSize" << " " << vpGS[k] << " "
														 << "dblXGridMin 500 dblYGridMin 194 "
														 << "R0_Total " << vpR0[l] << " "
														 << "intMaxGeneration " << vpMG[m] << " "
														 << "EndTime " << vpMT[n] << " "
														 << "RealisationsPerParameterSet " << vpRP[o] << " "
														 << "R0_Proportion_Spatial " << vpPS[p] << " "
														 << "Relative_Transmit_Prodromal " << vpTP[q] << " "
 														 << "Relative_Transmit_Neighbour_vs_House_Prodromal " << vpHP[r] << " "
 														 << "Relative_Transmit_Neighbour_vs_House_Symptomatic " << vpHS[s] << " "
														 << "\n";
												}
											}
										}
									}
								}
							}
						}					
					}
				}
			}
		}
	}
	ofs2 << "pause\n";
	ofs2.close();
	}
	if (arg=="ANALYSIS") {
	batchname = srjRES+"\\"+"analysis.txt";
	cerr << "Analysing results accross multiple files.  Writing output to " << batchname << "\n";
	vpPD.resize(0);
	vpPDn.resize(0);
	vpSD.resize(0);
	vpWS.resize(0);
	vpGS.resize(0);
	vpR0.resize(0);
	vpMG.resize(0);
	vpMT.resize(0);
	vpRP.resize(0);
	vpPS.resize(0);
	vpTP.resize(0);
	vpHP.resize(0);
	vpHS.resize(0);

	vpPD.push_back(srjDATA+"\\EWPopDense.txt "+srjDATA +"\\EWWorkDense.txt ");
	vpPD.push_back(srjDATA+"\\FlatDense.txt "+srjDATA +"\\FlatDense.txt ");
	vpPDn.push_back("RD");
	vpPDn.push_back("FD");
	vpSD.push_back("3984029522");
	// vpSD.push_back("4958948958");
	// vpSD.push_back("2309430938");
	// vpSD.push_back("5834729388");
	// vpSD.push_back("1928499182");
	// vpSD.push_back("0594900321");
	vpWS.push_back("5");
	// vpWS.push_back("10");
	// vpWS.push_back("15");
	vpWS.push_back("20");
	vpGS.push_back("6.392809287");
	vpGS.push_back("14.29475613");
	vpGS.push_back("20.21583799");
	vpGS.push_back("45.20398798");
	vpGS.push_back("63.92809287");
	vpR0.push_back("4");
	vpR0.push_back("6");
	vpR0.push_back("8");
	vpMG.push_back("4");
	vpMT.push_back("100");
	vpRP.push_back("1000");
	vpPS.push_back("1");
	vpPS.push_back("0.1");
	vpTP.push_back("0.2");
	vpHP.push_back("0.25");
	vpHS.push_back("0.125");

	pair<double,double> tmppair;

	ofs2.open(batchname.c_str());
	ofs2 << "PD\tSD\tWS\tGS\tR0\tMG\tMT\tRP\tPS\tTP\tHP\tHS\taveSPX\n";
	if (ofs2.fail()) {cerr << "Error opening output file " << batchname << "\n";exit(1);};
	for (int h=0;h<vpPD.size();++h) {
		for (int i=0;i<vpGS.size();++i) {
			for (int j=0;j<vpWS.size();++j) {
				for (int k=0;k<vpSD.size();++k) {
					for (int l=0;l<vpR0.size();++l) {
						for (int m=0;m<vpMG.size();++m) {
							for (int n=0;n<vpMT.size();++n) {
								for (int o=0;o<vpRP.size();++o) {
									for (int p=0;p<vpPS.size();++p) {
										for (int q=0;q<vpTP.size();++q) {
											for (int r=0;r<vpHP.size();++r) {
												for (int s=0;s<vpHS.size();++s) {
													ossf.str("");
													ossf << vpPDn[h] << "_GS_" << vpGS[i] << "_WS_" << vpWS[j] << "_SD_" << vpSD[k] << ".hex";
													osso.str("");
													osso << ossf.str() << "_R0_" << vpR0[l] << "_MG_" << vpMG[m] << "_MT_" << vpMT[n] << "_RP_" << vpRP[o] << "_PS_" << vpPS[p] << "_TP_" << vpTP[q] << "_HP_" << vpHP[r] << "_HS_" << vpHS[s] << "_SPX.txt";
													cerr << "Starting " << osso.str() << "...";
													srtp.LoadFromFile(srjRES+"\\"+osso.str());
													tmppair = srtp.GetLastPercentiles(0.025,0.975);
													ofs2 << vpPDn[h] << "\t" << vpGS[i] << "\t" << vpWS[j] << "\t" << vpSD[k] << "\t" << vpR0[l] << "\t" << vpMG[m] << "\t" << vpMT[n] << "\t" << vpRP[o] << "\t" << vpPS[p] << "\t" << vpTP[q] << "\t" << vpHP[r] << "\t" << vpHS[s] << "\t";
													ofs2 << srtp.GetLastAverage() << "\t" << tmppair.first << "\t" << tmppair.second << "\n";
													cerr << "done.\n";
												}
											}
										}
									}
								}
							}
						}					
					}
				}
			}
		}
	}
	ofs2 << "pause\n";
	ofs2.close();
	}
	return 0;
}

void SR::TemporalData::LoadFromFile(string filename) {
	ifstream ifs;
	ifs.open(filename.c_str());
	string junk;
	double tmpdouble;
	if (ifs.fail()) {cerr<<"Problem opening data file in SR::TemporalData::LoadFromFile";exit(1);};
	ifs >> junk >> noRealisations >> junk >> noTimeSteps;
	if ((noRealisations+2)*(noTimeSteps+1)>rawData.size()) {cerr<<"Data file too large in SR::TemporalData::LoadFromFile";exit(1);};
	for (int i=0;i<noTimeSteps+1;++i) ifs >> junk;
	for (int i=0;i<noRealisations+2;++i) {
		ifs >> rawData[i*(noTimeSteps+1)+0] >> junk;
		for (int j=0;j<noTimeSteps;++j) {
			ifs >> rawData[i*(noTimeSteps+1)+j+1];
		}
	}
	ifs.close();
}

double SR::TemporalData::GetLastAverage() {
	return rawData[noTimeSteps];
};

pair<double,double> SR::TemporalData::GetLastPercentiles(double lp, double up) {
	pair<double,double> rtnpr;
	int lbindex,ubindex;
	vector<double> lastvals(noRealisations,-1);
	for (int i=0;i<noRealisations;++i) {
		lastvals[i] = rawData[i*(noTimeSteps+1)+noTimeSteps];
	}
	sort(lastvals.begin(),lastvals.end());
	lbindex = static_cast<int>(ceil(1.0*noRealisations*lp));
	ubindex = static_cast<int>(floor(1.0*noRealisations*up));
	if (lbindex >= noRealisations) {
		cerr<<"index too high\n";
		exit(1);
		lbindex=noRealisations-1;
	}
	rtnpr.first = lastvals[lbindex];
	rtnpr.second = lastvals[ubindex];
	return rtnpr;
};