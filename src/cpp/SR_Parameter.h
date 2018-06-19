#ifndef SR_INC_PARAMETER
#define SR_INC_PARAMETER

#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include<algorithm>
#include<iomanip>
#include<time.h>
// #include"nr.h"
// #include"SR_Utility.h"
#include"SR_Label.h"
// #include"SR_Utility.h"
#include"SR_StripComments.h"

using namespace std;

namespace SR {

	void srerror(const string etext);

	class ParameterDash {
	private:
		bool numeric;
		Label name,tag;
		double value,initialvalue;
	public:
		double Value() {return value;};
		string Tag() {return tag.Get();}
		double* Pointer() {return &value;};
		double* InitialPointer() {return &initialvalue;};
		ParameterDash(): name(),tag() {numeric = true;};
		void Name(string s) {name.Set(s);};
		void Set(double d) {value = d;};
		void Set(string s) {tag.Set(s);};
		void SetName(string s) {name.Set(s);};
		string GetName() {return name.Get();};
		bool IsNumeric() {return numeric;};
		void SetNumeric(bool b) {numeric=b;};
		void Assign(string s, bool n, double d, string t);
		void SetInitialValue() {initialvalue=value;};
		void SetInitialValue(double d) {initialvalue=d;};
		double GetInitialValue() {return initialvalue;};
		void RevertToInitialValue() {value=initialvalue;};
		friend ofstream& operator<<(ofstream& ofs, ParameterDash& p);
		friend ifstream& operator>>(ifstream& ifs, ParameterDash& p);
	};

	class ParameterSet {
	private:
		static const int max = 150;
		int n;
		vector<ParameterDash> p;
		bool locked;
		vector<ParameterDash>::iterator GetParameterPointer(string s);
		vector<ParameterDash>::iterator FindName(string s);
	public:
		int intSeed;
		ParameterSet():p(max){n=0;locked=false;intSeed=-1234;};
		void AddValue(string s, double d);
		void AddValue(string s, string t);
		double* GetPointer(string s);
		double* GetInitialPointer(string s);
		double GetValue(string s);
		int GetIntValue(string s) {return static_cast<int>(GetValue(s));};
		string GetTag(string s);
		void ReadParams(string s);
		void ReadParamsFromFile(string s);
		string WriteParams();
		void Lock();
		void RevertToInitialValues();
		inline void ChangeValue(string s, double d) {double *p=GetPointer(s);*p=d;};
		inline void ChangeInitialValue(string s, double d) {vector<ParameterDash>::iterator it=GetParameterPointer(s);it->SetInitialValue(d);};
		friend ofstream& operator<<(ofstream& ofs, ParameterSet& ps);
		friend ifstream& operator>>(ifstream& ifs, ParameterSet& ps);
	};

	class ParameterValueSet {
	private:
		int noChanges;
		string **paramlabels;
		double *paramvalues;
		int *setnumbers;
		int currentparamnumber;
		int maxsetnumber;
	public:
		ParameterValueSet() {srerror("No default constructor for ParameterValueSet");};
		ParameterValueSet(string filename);
		~ParameterValueSet();
		void UpdateParameterSet(SR::ParameterSet &p, int ps_num);
		int MaxSetNumber() {return maxsetnumber;};
	};

	class ParameterValueSetB {
	private:
		int noChanges;
		int noParams;
		string **paramlabels;
		double *paramvalues;
		int currentparamnumber;
	public:
		ParameterValueSetB() {srerror("No default constructor for ParameterValueSet");};
		ParameterValueSetB(string filename);
		~ParameterValueSetB();
		void UpdateParameterSet(SR::ParameterSet &p);
		int NoChanges() {return noChanges;};
	};

	ParameterSet CommandLineAndParameterFile(int argc, char **argv);

}

#endif
