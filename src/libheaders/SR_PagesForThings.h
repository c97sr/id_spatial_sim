#ifndef SR_INC_PAGESFORTHINGS
#define SR_INC_PAGESFORTHINGS

#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include<algorithm>
#include<iomanip>
#include<time.h>
#include"nr.h"
#include"SR_Utility.h"

/*

 A template class to set up pages memory for lists of pointers to things of type T.
 On reflection - probably doesn't need to be templated - but it was a good exercise
 none the less.  The SR_PAE_PAGING needs to be set as a command line option to force
 the program to use greater than 3GB of memory.

 You need the windows core sdk installed to use "windows.h".  Look up some of the AWE
 commands at Microsoft.com for more info.

 Don't forget tne /D_WIN32_WINNT=0x0502 compile option to tell the compiler you're using
 Windows 2000 Advanced Server or better.

 Copy-paste command line options:
 /D_WIN32_WINNT=0x0502 /D SR_PAE_PAGING

 Page with AWE example

 http://msdn.microsoft.com/library/default.asp?url=/library/en-us/memory/base/awe_example.asp
 http://www.winnetmag.com/Article/ArticleID/7290/7290.html

 S Riley 20 July 2004

*/

#ifdef SR_PAE_PAGING
#include<windows.h>
BOOL LoggedSetLockPagesPrivilege ( HANDLE hProcess, BOOL bEnable);
#endif

using namespace std;

namespace SR {

	// Paged pointer
	template<class T> class PagedThingPointer {
		int intPageIndex,ptThing;
	public:
		inline int GetPageIndex();
		inline int GetPtThing();
		inline void SetPageIndex(int i);
		inline void SetPtThing(int i);
	};

	// Paging for nodes
	template<class T> class PagesForThings {
	private:
		unsigned int intCurrentBlockUsed;
		static const int intInternalPageSize = 4096;
		T** ptrT;
		unsigned int intNoBlocks;
		int intBlockSize;
		int currentPage;
		vector<int> vecNumberOfThingPointers;
		T** FirstOfPage(int i);
		vector<T**> ptFirstThing;
#ifdef SR_PAE_PAGING
		int intBlockSizeInInternalPages;
		T** ptrStartWindow;
		BOOL bResult;                   // generic Boolean value
		ULONG_PTR NumberOfPages;        // number of pages to request
		ULONG_PTR NumberOfPagesInitial; // initial number of pages requested
		ULONG_PTR *aPFNs;               // page info; holds opaque data
		PVOID lpMemReserved;            // AWE window
		SYSTEM_INFO sSysInfo;           // useful system information
		int PFNArraySize;               // memory to request for PFN array
#endif
	public:
		PagesForThings<T>(int size, int no);
		PagesForThings<T>(int size, int no, ifstream& ifs, T* v);
		~PagesForThings();
		T** ReturnThing(PagedThingPointer<T> p);
		PagedThingPointer<T> InsertThing(int noNewThings);
		void WriteToBinaryFile(ofstream& ofs, T* v);
	};

}

template <class T> SR::PagesForThings<T>::PagesForThings(int size, int no) : vecNumberOfThingPointers(no,0), ptFirstThing(no) {
	intCurrentBlockUsed = 0;
	intNoBlocks=no;
	intBlockSize=size;
#ifndef SR_PAE_PAGING
	for (unsigned int i=0;i<intNoBlocks;++i) {
		ptrT=static_cast<T**>(calloc(intBlockSize,sizeof(T**)));
		if(!ptrT) SR::srerror("Error with calloc.\n");
		ptFirstThing[i]=ptrT;
	}
#else
	GetSystemInfo(&sSysInfo);  // fill the system information structure
	if (sSysInfo.dwPageSize != intInternalPageSize) SR::srerror("This version of the code is compiled for 32-bit PAE with specific internal page size.");
	if (intBlockSize % intInternalPageSize != 0 ) SR::srerror("Block size must be a multiple of internal page size.");
	intBlockSizeInInternalPages = intBlockSize * sizeof(T*) / intInternalPageSize;
	NumberOfPages = intNoBlocks*intBlockSizeInInternalPages;
	PFNArraySize = NumberOfPages * sizeof (ULONG_PTR);
	aPFNs = (ULONG_PTR *) HeapAlloc (GetProcessHeap (), 0, PFNArraySize);
	if (aPFNs == NULL) SR::srerror("Failed to allocate page pointers on heap.");
	if (!LoggedSetLockPagesPrivilege(GetCurrentProcess(),TRUE)) SR::srerror("Couldn't set the correct privilage for this user.");
	NumberOfPagesInitial = NumberOfPages;
	bResult = AllocateUserPhysicalPages(GetCurrentProcess(),&NumberOfPages,aPFNs);
	if(bResult!=TRUE || NumberOfPagesInitial != NumberOfPages) SR::srerror("Couldn't allocate physical pages.");
	lpMemReserved = VirtualAlloc( NULL,intBlockSize*sizeof(T*),MEM_RESERVE | MEM_PHYSICAL,PAGE_READWRITE);
	if(lpMemReserved == NULL) SR::srerror("Cannot reserve virtual memory in the process address space.");
	bResult = MapUserPhysicalPages(lpMemReserved,intBlockSizeInInternalPages,aPFNs);
	if (!bResult) {
		SR::srerror("Problem mapping physical pages");
	}
	ptrStartWindow = static_cast<T**>(lpMemReserved);
#endif
	currentPage=0;
};

template <class T> T** SR::PagesForThings<T>::ReturnThing(SR::PagedThingPointer<T> p) {
#ifndef SR_PAE_PAGING
	if (currentPage != p.GetPageIndex()) currentPage=p.GetPageIndex();
	return (ptFirstThing[currentPage]+p.GetPtThing());
#else
	if (currentPage != p.GetPageIndex()) {
		currentPage = p.GetPageIndex();
		bResult = MapUserPhysicalPages(lpMemReserved,intBlockSizeInInternalPages,aPFNs+currentPage*intBlockSizeInInternalPages);
		if (bResult != TRUE) {
			SR::srerror("Problem mapping physical pages");
		}
	}
	return ptrStartWindow + p.GetPtThing();
#endif
};

template <class T> SR::PagedThingPointer<T> SR::PagesForThings<T>::InsertThing(int noNewThings) {
	static unsigned int i;
	// problem may be here - sizeof needed?
	if (noNewThings > intBlockSize) SR::srerror("Block size too small.\n");
	PagedThingPointer<T> rtnPt;
	// cout << sizeof(rtnPt) << "\n";
	if (intCurrentBlockUsed < 100) i=0;
	else i = intCurrentBlockUsed-100;
	while (vecNumberOfThingPointers[i]+noNewThings>=intBlockSize) i++; // devide sizeof pointer?
	if (i==intNoBlocks) SR::srerror("Not enough higher memory reserved.\n");
	if (i>intCurrentBlockUsed) intCurrentBlockUsed=i;
	rtnPt.SetPageIndex(i);
	rtnPt.SetPtThing(vecNumberOfThingPointers[i]);
	vecNumberOfThingPointers[i] += noNewThings;
	return rtnPt;
};

template <class T> void SR::PagesForThings<T>::WriteToBinaryFile(ofstream& ofs, T* v) {
	T **ptptStart, **ptptEnd;
	int tmpint;
	cerr << "XXX" << intCurrentBlockUsed << " " << intNoBlocks << " " << intBlockSize << "\n";
	SR::BinWrite(ofs,intCurrentBlockUsed);
	SR::BinWrite(ofs,intNoBlocks);
	SR::BinWrite(ofs,intBlockSize);
	for (unsigned int i=0;i<intNoBlocks;++i) SR::BinWrite(ofs,vecNumberOfThingPointers[i]);
	for (unsigned int i=0;i<intNoBlocks;++i) {
		// This has to be the problem... no
		ptptStart = FirstOfPage(i);
		ptptEnd = ptptStart + vecNumberOfThingPointers[i];
		while (ptptStart != ptptEnd) {
			tmpint = *ptptStart-v;
			SR::BinWrite(ofs,tmpint);
			ptptStart++;
		}
	}
};

template <class T> T** SR::PagesForThings<T>::FirstOfPage(int index) {
	T** rtnval;
	if (index < 0 || index >= static_cast<int>(ptFirstThing.size()))  SR::srerror("Range problems in template<class T> T** SR::PagesForThings<T>::FirstOfPage(int index)");
#ifndef SR_PAE_PAGING
	rtnval = ptFirstThing[index];
	return rtnval;
#else
	if (index != currentPage) {
		currentPage = index;
		bResult = MapUserPhysicalPages(lpMemReserved,intBlockSizeInInternalPages,aPFNs+currentPage*intBlockSizeInInternalPages);
		if (bResult != TRUE) {
			SR::srerror("Problem mapping physical pages");
		}
	}
	return static_cast<T**>(ptrStartWindow);
#endif
};

template <class T> SR::PagesForThings<T>::PagesForThings(int size, int no, ifstream& ifs, T* v) :
vecNumberOfThingPointers(no,0), ptFirstThing(no) {
	T **ptptStart,**ptptEnd;
	int tmpint;
	intNoBlocks=no;
	intBlockSize=size;
	tmpint = SR::BinRead<int>(ifs);
	intCurrentBlockUsed = tmpint;
	cerr << tmpint << "\n";
	tmpint = SR::BinRead<int>(ifs);
	cerr << tmpint << "\n";
	if (tmpint != static_cast<int>(intNoBlocks)) SR::srerror("Network inconsistency in template<class T> SR::PagesForThings<T>::PagesForThings<T>(int size, int no, ifstream& ifs, vector<T>::iterator v)");
	tmpint = SR::BinRead<int>(ifs);
	cerr << tmpint << "\n";
	if (tmpint != intBlockSize) SR::srerror("Network inconsistency in template<class T> SR::PagesForThings<T>::PagesForThings<T>(int size, int no, ifstream& ifs, vector<T>::iterator v)");
	for (unsigned int i=0;i<intNoBlocks;++i) vecNumberOfThingPointers[i]=SR::BinRead<int>(ifs);
	cerr << "+Allocating memory in PagesForThings constructor (using binary file)...";
#ifndef SR_PAE_PAGING
	for (unsigned int i=0;i<intNoBlocks;++i) {
		ptrT=static_cast<T**>(calloc(intBlockSize,sizeof(T**)));
		if(!ptrT) SR::srerror("Error with calloc.\n");
		ptFirstThing[i]=ptrT;
	}
#else
	GetSystemInfo(&sSysInfo);  // fill the system information structure
	if (sSysInfo.dwPageSize != intInternalPageSize) SR::srerror("This version of the code is compiled for 32-bit PAE with specific internal page size.");
	if (intBlockSize % intInternalPageSize != 0 ) SR::srerror("Block size must be a multiple of internal page size.");
	intBlockSizeInInternalPages = intBlockSize * sizeof(T*) / intInternalPageSize;
	NumberOfPages = intNoBlocks*intBlockSizeInInternalPages;
	PFNArraySize = NumberOfPages * sizeof (ULONG_PTR);
	aPFNs = (ULONG_PTR *) HeapAlloc (GetProcessHeap (), 0, PFNArraySize);
	if (aPFNs == NULL) SR::srerror("Failed to allocate page pointers on heap.");
	if (!LoggedSetLockPagesPrivilege(GetCurrentProcess(),TRUE)) SR::srerror("Couldn't set the correct privilage for this user.");
	NumberOfPagesInitial = NumberOfPages;
	bResult = AllocateUserPhysicalPages(GetCurrentProcess(),&NumberOfPages,aPFNs);
	if(bResult!=TRUE || NumberOfPagesInitial != NumberOfPages) SR::srerror("Couldn't allocate physical pages.");
	lpMemReserved = VirtualAlloc( NULL,intBlockSize*sizeof(T*),MEM_RESERVE | MEM_PHYSICAL,PAGE_READWRITE);
	if(lpMemReserved == NULL) SR::srerror("Cannot reserve virtual memory in the process address space.");
	bResult = MapUserPhysicalPages(lpMemReserved,intBlockSizeInInternalPages,aPFNs);
	if (bResult != TRUE) {
		SR::srerror("Problem mapping physical pages");
	}
	ptrStartWindow = static_cast<T**>(lpMemReserved);
#endif
	cerr << "done.\n";
	cerr << "+Reading in network from binary file...\n";
	currentPage=0;
	for (unsigned int i=0;i<intNoBlocks;++i) {
		ptptStart = FirstOfPage(i);
		ptptEnd = ptptStart + vecNumberOfThingPointers[i];
		while (ptptStart != ptptEnd) {
			tmpint = SR::BinRead<int>(ifs);
			*ptptStart = v + tmpint;
			ptptStart++;
		}
		if (i%100==0) cerr << "+Completed block: " << i+1 << " of " << intNoBlocks << "                        \r";
	}
	cerr << "\n+...done.\n";
};

template <class T> SR::PagesForThings<T>::~PagesForThings() {
#ifndef SR_PAE_PAGING
	for (unsigned int i=0;i<intNoBlocks;++i) free(ptFirstThing[i]);
#else
	bResult = MapUserPhysicalPages(lpMemReserved,intBlockSizeInInternalPages,NULL);
	if (!bResult) SR::srerror("Failed to unmap pages in process address space");
	bResult = FreeUserPhysicalPages(GetCurrentProcess(),&NumberOfPages,aPFNs );
	if (!bResult) SR::srerror("Failed to unfree physical pages");
	bResult = VirtualFree( lpMemReserved,0,MEM_RELEASE );
	bResult = HeapFree(GetProcessHeap(), 0, aPFNs);
	if (!bResult) SR::srerror("Call to heap free has failed.");
#endif
};

template <class T> int SR::PagedThingPointer<T>::GetPageIndex() {
	return intPageIndex;
};

template <class T> int SR::PagedThingPointer<T>::GetPtThing() {
	return ptThing;
};

template <class T> void SR::PagedThingPointer<T>::SetPageIndex(int i) {
	intPageIndex=i;
};

template <class T> void SR::PagedThingPointer<T>::SetPtThing(int i) {
	ptThing=i;
};

#endif
