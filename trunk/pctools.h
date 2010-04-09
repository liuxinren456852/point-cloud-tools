// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the PCTOOLS_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// PCTOOLS_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef PCTOOLS_EXPORTS
#define PCTOOLS_API __declspec(dllexport)
#else
#define PCTOOLS_API __declspec(dllimport)
#endif

// This class is exported from the pctools.dll
//extern "C" class PCTOOLS_API Cpctools {
//public:
//	Cpctools(void);
	// TODO: add your methods here.
//};
/*
extern "C" {
	PCTOOLS_API int npctools;
	PCTOOLS_API int fnpctools(void);
	PCTOOLS_API int importpts(const _TCHAR* from_ascii, const _TCHAR* save);
	struct PCDescr;
	PCTOOLS_API int describe(const _TCHAR* file, PCDescr &descr);
	PCTOOLS_API int sortintoxygrid(const _TCHAR* file, double xres, double yres);
}*/
#include "tchar.h"
extern "C" PCTOOLS_API int NUM_THREADS;

// extern "C" PCTOOLS_API const unsigned NPOINTS			= 1<<0;
// extern "C" PCTOOLS_API const unsigned ZMEAN			= 1<<1;
// extern "C" PCTOOLS_API const unsigned ZMEAN_DETR		= 1<<2;
// extern "C" PCTOOLS_API const unsigned ZMIN			= 1<<3;
// extern "C" PCTOOLS_API const unsigned ZMAX			= 1<<4;
// extern "C" PCTOOLS_API const unsigned STDEV			= 1<<5;
// extern "C" PCTOOLS_API const unsigned STDEV_DETR		= 1<<6;
// extern "C" PCTOOLS_API const unsigned SKEW			= 1<<7;
// extern "C" PCTOOLS_API const unsigned SKEW_DETR		= 1<<8;
// extern "C" PCTOOLS_API const unsigned KURT			= 1<<9;
// extern "C" PCTOOLS_API const unsigned KURT_DETR		= 1<<10;

extern "C" PCTOOLS_API int set_num_threads(int num_threads);
extern "C" PCTOOLS_API int showconsole();
#include "define.h"
/*
Python:
class PCDESCR(Structure):
	_fields_ = [("npoints", c_ulonglong),
	("floatnbits", c_int),
	("pmin", c_double * 3),
	("pmax", c_double * 3),
	("pcom", c_double * 3),
	("p0", c_double * 3),
	("res", c_double * 3)
	]


descr = PCDESCR()

pctools.importpts(pts, unsorted_mmf, byref(descr))
*/
/*
struct PCDESCR {
index npoints;
int floatnbits; //64 or 32
double3 pmin, pmax, pcom, p0, res;
};
extern "C" PCTOOLS_API int importpts(const _TCHAR* from_ascii, const _TCHAR* save, PCDESCR &descr);

struct Rational {
	int num;
	int den;
};*/
extern "C" PCTOOLS_API int importpts(const _TCHAR* from_ascii, const _TCHAR* save);

//extern "C" struct MMFini;
//extern "C" PCTOOLS_API int readmmfini(const _TCHAR* file, MMFini &ini);
//extern "C" PCTOOLS_API int writemmfini(const _TCHAR* file, const MMFini &ini);
extern "C" PCTOOLS_API int sortintoxygrid(const _TCHAR* unsortedpcmmf_filename, uint32 xres_nom, uint32 yres_nom, uint32 denominator, double normal[3]);
//extern "C" PCTOOLS_API int LSPCA(const _TCHAR* pcmmf_filename);
extern "C" PCTOOLS_API int changegridresolution(const _TCHAR* pcmmf_filename, uint32 xres_nom, uint32 yres_nom, uint32 denom);

//extern "C" PCTOOLS_API int stats(const _TCHAR* pcmmf_filename);
extern "C" PCTOOLS_API int stats_detrended(const _TCHAR* pcmmf_filename, const double times_zmin, const double times_zmean, const double times_stdev);
extern "C" PCTOOLS_API int export_stats(const _TCHAR* pcmmf_filename, const uint64 min_npoints_per_cell);

extern "C" PCTOOLS_API int exportply(const _TCHAR* mmf_name, const _TCHAR* ply_name);

extern "C" PCTOOLS_API int rasterizeongrid(
					const _TCHAR* pc_mmf_name,
					const _TCHAR* grid_mmf_name
					);

extern "C" PCTOOLS_API int rasterizeonsubgrid(
					   const _TCHAR* pc_mmf_name,
					   const _TCHAR* grid_mmf_name,
					   const unsigned int nsubd,
					   const _TCHAR* subgrid_mmf_name,
					   const _TCHAR* detrended_subgrid_mmf_name
					   );


//extern "C" struct MMFini {
//	uint64 count, nx, ny;
//	uint32 floatnbits, subd;
//	double xmin, xmax, ymin, ymax, zmin, zmax, xres0, yres0;
//};
//const MMFini ini_null = {
//	0, 0, 0, 0, 0,
//	NaN, NaN, NaN, NaN, NaN, NaN, 0., 0., //8
//};

