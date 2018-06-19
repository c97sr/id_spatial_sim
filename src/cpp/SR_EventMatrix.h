#ifndef SR_INC_EVENTMATRIX
#define SR_INC_EVENTMATRIX

#include<numeric>
#include"SR_Utility.h"
#include"SR_GatherRunInformation.h"
#include"SR_GridHex.h"

/*

S Riley 

*/

using namespace std;

namespace SR {
	// Event matrix
	class EventMatrix {		
	private:
		int intMaxDelay,intMaxEvents,intCurrentMatrixTime,intLevelConstraint,intTotalEvents;
		SR::UNTIMEDEVENT constrainedEv;
		vector<SR::UntimedEvent> matStacksLinear;
		vector<int> vecSizeOfStacks;
		vector<int> vecNumberOfConstrainedEvents;
	public:
		EventMatrix() {SR::srerror("No default constructor for EventMatrix\n");};
		EventMatrix(int md, int me,UNTIMEDEVENT cev,int levcon);
		void AddEvent(UntimedEvent ue, int delay);
		void ApplyEvents(SR::ParameterSet& p);
		void ApplyEvents(SR::ParameterSet& p, GatherRunInformation& ri);
		int TotalEventsPending();
		int EventsPendingAtThisTimeStep();
		void ClearAllEvents();
		string StringOutputTheStack();
		void RemoveDuplicatesAtCurrentTime();
		bool CurrentTimeHasEvents();
		UntimedEvent PopTop();
		int GetMaxDelay() {return intMaxDelay;};
		void ExportEvents(SR::EventMatrix& tmpEm);
		void RemoveEventAtRandom(SR::ParameterSet& p);
		void SaveStackToFile(string filename, SR::GridHex& gh, int (*f1)(SR::UNTIMEDEVENT ue));
		void LoadPendingEventsFromFile(string filename, SR::GridHex& gh,SR::UNTIMEDEVENT (*f1)(int));
	};
}

#endif
