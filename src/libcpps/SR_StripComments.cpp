#include"SR_StripComments.h"

string SR::StripComments(string s) {
	ostringstream oss;
	string::iterator p=s.begin();
	bool copyon=true;
	while (p != s.end()) {
		if (copyon && *p=='[') copyon=false;
		if (copyon) oss << *p;
		if (!copyon && *p==']') copyon=true;
		p++;
	}
	return oss.str();
};
