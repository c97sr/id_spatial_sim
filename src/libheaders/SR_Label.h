#ifndef SR_INC_LABEL
#define SR_INC_LABEL

#include<iomanip>
#include<vector>
#include"nr.h"
// #include"SR_Utility.h"

using namespace std;

namespace SR {

	void srerror(const string etext);

	class Label {
	private:
		static const int noChars = 512;
		int size;
		vector<char> data;
	public:
		Label():data(noChars){size=0;};
		void Set(string s);
		Label(string s) : data(noChars) {Set(s);};
		string Get();
		int GetSize() {return size;};
		bool IsEqual(string s);
		friend ofstream& operator<<(ofstream& ofs, Label& lab);
		friend ifstream& operator>>(ifstream& ifs, Label& lab);
	};
	ofstream& operator<<(ofstream& ofs, SR::Label& lab);
	ifstream& operator>>(ifstream& ifs, SR::Label& lab);
}

#endif
