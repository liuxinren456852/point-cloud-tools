// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// Modify the following defines if you have to target a platform prior to the ones specified below.
// Refer to MSDN for the latest info on corresponding values for different platforms.
#ifndef WINVER				// Allow use of features specific to Windows XP or later.
#define WINVER 0x0501		// Change this to the appropriate value to target other versions of Windows.
#endif

#ifndef _WIN32_WINNT		// Allow use of features specific to Windows XP or later.                   
#define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other versions of Windows.
#endif						

#ifndef _WIN32_WINDOWS		// Allow use of features specific to Windows 98 or later.
#define _WIN32_WINDOWS 0x0410 // Change this to the appropriate value to target Windows Me or later.
#endif

#ifndef _WIN32_IE			// Allow use of features specific to IE 6.0 or later.
#define _WIN32_IE 0x0600	// Change this to the appropriate value to target other versions of IE.
#endif

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:
//#include <windows.h>



// TODO: reference additional headers your program requires here
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS // This is not enough sometimes, 
#pragma warning( disable : 4996)

#include <io.h>
//#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <tchar.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstddef>
#include <cstdio>    //std::remove
#include <iostream>


#include <omp.h>


#include <ios>
#define FUSION_MAX_VECTOR_SIZE 15
//#define FUSION_MAX_MAP_SIZE 12
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/typeof/typeof.hpp>

#define BOOST_DATE_TIME_NO_LIB
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <boost/format.hpp>
#include <boost/format.hpp>

#include <boost/range.hpp>

#include <boost/assert.hpp>
//#include <Eigen/Core>

//#include "define.h"
using namespace std;

