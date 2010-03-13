// pctools.cpp : Defines the entry point for the DLL application.
//
#include "stdafx.h"
#include <windows.h>
#include "pctools.h"


#ifdef _MANAGED
#pragma managed(push, off)
#endif

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
    return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif

// This is an example of an exported variable
PCTOOLS_API int NUM_THREADS=8;

// This is an example of an exported function.
PCTOOLS_API int set_num_threads(int num_threads)
{
	NUM_THREADS=num_threads;
	return 0;
}

// This is the constructor of a class that has been exported.
// see pctools.h for the class definition
//Cpctools::Cpctools()
//{
//	return;
//}


int showconsole() {
	BOOL bres = AllocConsole();// allocate a console for this app
//	if(!bres && 5 == GetLastError()) {
	if(!bres) {
//		cout << "we already have a console?";
		FreeConsole();
		/*BOOL bres = */AllocConsole();// allocate a console for this app
//		return 0; // we already have a console
	}
	// redirect unbuffered STDOUT to the console
	intptr_t lStdHandle = (intptr_t)GetStdHandle(STD_OUTPUT_HANDLE);
	int hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
	if(hConHandle == -1)
		return 1;
	FILE *fp = _fdopen( hConHandle, "w" );
	*stdout = *fp;
	setvbuf( stdout, NULL, _IONBF, 0 );
	 // redirect unbuffered STDERR to the console
	lStdHandle = (intptr_t)GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
	fp = _fdopen( hConHandle, "w" );
	*stderr = *fp;
	setvbuf( stderr, NULL, _IONBF, 0 );
	// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog point to console as well
	ios::sync_with_stdio();
	//cout << "cout "; cerr << "cerr "; clog << "clog "; printf("printf\n");
	return 0;
}