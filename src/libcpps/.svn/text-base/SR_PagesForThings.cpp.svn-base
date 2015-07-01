#include"SR_PagesForThings.h"

/*

S Riley 20th July 2004

*/

#ifdef SR_PAE_PAGING
BOOL LoggedSetLockPagesPrivilege ( HANDLE hProcess, BOOL bEnable) {
	struct { DWORD Count;LUID_AND_ATTRIBUTES Privilege [1];} Info;
	HANDLE Token;
	BOOL Result;

	// Open the token.
	Result = OpenProcessToken ( hProcess,
		TOKEN_ADJUST_PRIVILEGES,
		& Token);

	if( Result != TRUE ) 
	{
		cerr << "Cannot open process token.\n";
		return FALSE;
	}

	// Enable or disable?
	Info.Count = 1;
	if( bEnable ) 
	{
		Info.Privilege[0].Attributes = SE_PRIVILEGE_ENABLED;
	} 
	else 
	{
		Info.Privilege[0].Attributes = 0;
	}

	// Get the LUID.
	Result = LookupPrivilegeValue ( NULL,
		SE_LOCK_MEMORY_NAME,
		&(Info.Privilege[0].Luid));

	if( Result != TRUE ) 
	{
		cerr << "Cannot get privilege value for " << SE_LOCK_MEMORY_NAME << ".\n";
		return FALSE;
	}

	// Adjust the privilege.
	Result = AdjustTokenPrivileges ( Token, FALSE,
		(PTOKEN_PRIVILEGES) &Info,
		0, NULL, NULL);

	// Check the result.
	if( Result != TRUE ) 
	{
		cerr << "Cannot adjust token privileges, error " << GetLastError() << ".\n";
		return FALSE;
	} 
	else 
	{
		if( GetLastError() != ERROR_SUCCESS ) 
		{
			cerr << "Cannot enable SE_LOCK_MEMORY_NAME privilege, please check the local policy.\n";
			return FALSE;
		}
	}
	CloseHandle( Token );
	return TRUE;
};

#endif
