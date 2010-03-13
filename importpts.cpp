#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmfini.h"

//int importpts(const _TCHAR* from_ascii, const _TCHAR* mmf, PCDESCR & descr)
int importpts(const _TCHAR* from_ascii, const _TCHAR* mmf)
{
	time_t start = time(NULL);
	uint64 nl = std::numeric_limits<uint64>::max();
	FILE *xyz_text = _tfopen(from_ascii, _T("rtS"));
	FILE *xyz_binary = _tfopen(mmf, _T("wbS"));
	if(xyz_text == NULL || xyz_binary == NULL) {
      cerr << "Problem opening file\n";
	  return 1;
	}
	double3 pmin = {DBL_MAX,DBL_MAX,DBL_MAX}, pmax={-DBL_MAX,-DBL_MAX,-DBL_MAX}, pcom = {0,0,0};

	uint64 l = 0,
		npoints = 0; // number of valid points

	char line[512];
	fgets(line, sizeof(line), xyz_text);// Always skipping the first line
	printf("skippping line #%llu: %s", l, line);
	l++;

	while(l < nl) {
		char*str = fgets(line, sizeof(line), xyz_text);
		if(str == NULL) //ferror(xyz_text) || feof(xyz_text))
			break; //assuming EOF
		if(*str == 0) {//ferror(xyz_text) || feof(xyz_text))
			cerr << "Corrupted file, binary 0 at the line # " << l << endl;
			return 1;
		}
		l++;
		double3 p;
		int nfields = sscanf(line, "%lf %lf %lf", &p.x, &p.y, &p.z);
		if(nfields != 3) {
			fprintf(stderr, "skippping line #%llu: %s", l, line);
			continue;
		}

		auto update_limits_along = [&p, &pmin, &pmax](int coord) {
			auto& x = as_array(p)[coord]; auto& xmin = as_array(pmin)[coord]; auto& xmax = as_array(pmax)[coord];
			if(x < xmin) xmin = x;
			if(x > xmax) xmax = x;
		};
		update_limits_along(X);
		update_limits_along(Y);
		update_limits_along(Z);

		pcom+=p;
		npoints++;
		fwrite(as_array(p), sizeof(p), 1, xyz_binary);
		static_assert(sizeof(p)==3*8, "wrong size");
		if(l%(1000000)==0)
			printf("\r%llu", l);
	}
	fclose(xyz_binary);
	fclose(xyz_text);
	pcom /= npoints;

	// min and max must encompass all the points, even for the truncated floating point precision
// 	pmin.x = min<double>(pmin.x, float(pmin.x));
// 	pmin.y = min<double>(pmin.y, float(pmin.y));
// 	pmin.z = min<double>(pmin.z, float(pmin.z));
// 	pmax.x = max<double>(pmax.x, float(pmax.x));
// 	pmax.y = max<double>(pmax.y, float(pmax.y));
// 	pmax.z = max<double>(pmax.z, float(pmax.z));

	mmfini::pc ini;
	uint32 floatnbits = sizeof(double)*8;
	const double3 p0 = {0,0,0};
#define ASSIGN(field) bf::at_key<mmfini::field>(ini) = field
	ASSIGN(npoints);
	ASSIGN(floatnbits);
	ASSIGN(pmin);
	ASSIGN(pmax);
	ASSIGN(pcom);
	ASSIGN(p0);
#undef ASSIGN

	const wstring filename_ext(mmf);
	const wstring::size_type beginext = filename_ext.rfind('.');
	const wstring filename(filename_ext, 0, beginext);
//	const wstring ext(filename_ext, beginext);

	if(!mmfini::save(ini, mmfini::filename(filename).c_str())) {
		cerr << "Cannot save description file";
		return 1;
	}


//	PCDESCR descr1 = {npoints, 64, pmin, pmax, pcom, p0, {0,0,0}};
//	descr = descr1;

	clog << "\r" << l << "lines (" << npoints << " points) in " << time(NULL) - start << "sec\n";
	return 0;
}