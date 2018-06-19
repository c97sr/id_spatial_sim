#ifndef SR_INC_GATHERRUNINFORMATION
#define SR_INC_GATHERRUNINFORMATION

#include<iostream>
// #include"nr.h"
#include"SR_Utility.h"

/*
 
[Enter details]

S Riley [Enter details]

*/

using namespace std;

namespace SR {
	class GatherRunInformation {
	public:
		virtual void RegisterEventAfterApplication(vector<UntimedEvent>::iterator ptEv, int cb)=0;
		virtual string ConvertEventPointerToString(UNTIMEDEVENT ue)=0;
	};
}

#endif
