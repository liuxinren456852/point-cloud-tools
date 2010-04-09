#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmfini.h"

/*
struct file_closer {
	void operator()(FILE* f) const { if (f) std::fclose(f); }
};

typedef unique_ptr<FILE,file_closer> file_handle;
*/


struct scope_guard_file_closer {
	FILE* f;
	scope_guard_file_closer(FILE* fin) : f(fin) {}
	~scope_guard_file_closer() { if (f) fclose(f); }
};


/*
inline void KahanAdd(double3 &sum, double3 p) {
	T reduceCPU(T *data, int size)
	{
		double3 sum = data[0];
		double3 c = {0,0,0};              
		for (int i = 1; i < size; i++)
		{
			T y = data[i] - c;  
			T t = sum + y;      
			c = (t - sum) - y;  
			sum = t;            
		}
		return sum;
}
*/
//int importpts(const _TCHAR* from_ascii, const _TCHAR* mmf, PCDESCR & descr)
int importpts(const _TCHAR* from_ascii, const _TCHAR* mmf)
{
	time_t start = time(NULL);
	uint64 nl = std::numeric_limits<uint64>::max();
	FILE *xyz_text = _tfopen(from_ascii, _T("rtS"));
	scope_guard_file_closer sgfc1(xyz_text);
	FILE *xyz_binary = _tfopen(mmf, _T("wbS"));
	scope_guard_file_closer sgfc2(xyz_binary);
	if(xyz_text == NULL || xyz_binary == NULL) {
      cerr << "Problem opening file\n";
	  return 1;
	}

	uint64 l = 0,
		npoints = 0; // number of valid points

	char line[512];
	fgets(line, sizeof(line), xyz_text);// Always skipping the first line
	printf("skippping line #%llu: %s", l, line);
	l++;


	double3 p0 = {0,0,0};
	fpos_t pos;
	if(fgetpos(xyz_text,&pos)) {
		cerr << "fgetpos error\n";
		return 1;
	}
	while(l < nl) {
		char*str = fgets(line, sizeof(line), xyz_text);
		if(str == NULL || *str == 0) {//ferror(xyz_text) || feof(xyz_text))
			cerr << "Empty or Corrupted file\n";
			return 1;
		}
		l++;
		double3 p;
		int nfields = sscanf(line, "%lf %lf %lf", &p.x, &p.y, &p.z);
		if(nfields != 3) {
			fprintf(stderr, "skippping line #%llu: %s", l, line);
		} else {
			p0 = p;
			p-=p0;
			fwrite(as_array(p), sizeof(p), 1, xyz_binary);
			break; // exit
		}
	}
	//double3 pmin = {DBL_MAX,DBL_MAX,DBL_MAX}, pmax={-DBL_MAX,-DBL_MAX,-DBL_MAX}, pcom = {0,0,0};
	npoints = 1;
	double3 pmin = {0,0,0}, pmax = {0,0,0};
	double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
	FPSUM(xx);
	FPSUM(xy);
	FPSUM(xz);
	FPSUM(yy);
	FPSUM(yz);
	FPSUM(zz);
	double3 pcom = {0,0,0};
	FPSUM(pcom);

//	if(!fsetpos(xyz_text,&pos)) {
//		cerr << "fsetpos error\n";
//		return 1;
//	}


	//auto &sum = pcom;
	//FPADDINIT(pcom);
	//KahanSum<decltype(pcom)> pcom_ks(pcom);
//	FPADDINIT(pcom);
	//FPSimpleAdd<decltype(pcom)> pcom_add(pcom);


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

		p-=p0;

		auto update_limits_along = [&p, &pmin, &pmax](int coord) {
			auto& x = as_array(p)[coord]; auto& xmin = as_array(pmin)[coord]; auto& xmax = as_array(pmax)[coord];
			if(x < xmin) xmin = x;
			if(x > xmax) xmax = x;
		};
		update_limits_along(X);
		update_limits_along(Y);
		update_limits_along(Z);

		const double x = p.x;
		const double y = p.y;
		const double z = p.z;
		xx_add(x*x);
		xy_add(x*y);
		xz_add(x*z);
		yy_add(y*y);
		yz_add(y*z);
		zz_add(z*z);
		pcom_add(p);

		fwrite(as_array(p), sizeof(p), 1, xyz_binary);
		npoints++;


		if(l%(1000000)==0)
			printf("\r%llu", l);
	}
	//fclose(xyz_binary);
	//fclose(xyz_text);
	pcom /= npoints;
	auto const x = pcom.x, y = pcom.y, z = pcom.z;
	xx -= x*x * npoints;
	xy -= x*y * npoints;
	xz -= x*z * npoints; 
	yy -= y*y * npoints;
	yz -= y*z * npoints;
	zz -= z*z * npoints;

	Wm4::Matrix3d const	A(
		xx,xy,xz,
		xy,yy,yz,
		xz,yz,zz);

	//clog << "A = \n" << A << endl;

	// Bounded-time eigensolver.
	Wm4::NoniterativeEigen3x3d const eigen(A);

// 	clog << "eigen values: \n";
// 	clog << eigen.GetEigenvalue(0) << " ";
// 	clog << eigen.GetEigenvalue(1) << " ";
// 	clog << eigen.GetEigenvalue(2) << "\n";
// 	clog << "eigen values: \n";
// 	clog << eigen.GetEigenvector(0) << "\n";
// 	clog << eigen.GetEigenvector(1) << "\n";
// 	clog << eigen.GetEigenvector(2) << "\n";

// 	Wm4::Vector3d const
// 		normal = eigen.GetEigenvector(0),
// 		ez(0,0,1),
// 		ey(0,1,0);
// 
// 	Wm4::Quaterniond q1,q2;
// 	q1.Align(normal, ez);
// 
// 	Wm4::Vector3d slope = q1.Rotate(-ez + normal*normal.Z());
// 	slope.Normalize();
// 	q2.Align(slope, -ey);
// 
// 	Wm4::Quaterniond const q = q2*q1;

// 	clog << "q = " << q << endl;
// 
// 	Wm4::Vector3d const
// 		normal_new = q.Rotate(normal),
// 		slope_new = q.Rotate(-ez + normal*normal.Z());


	// min and max must encompass all the points, even for the truncated floating point precision
// 	pmin.x = min<double>(pmin.x, float(pmin.x));
// 	pmin.y = min<double>(pmin.y, float(pmin.y));
// 	pmin.z = min<double>(pmin.z, float(pmin.z));
// 	pmax.x = max<double>(pmax.x, float(pmax.x));
// 	pmax.y = max<double>(pmax.y, float(pmax.y));
// 	pmax.z = max<double>(pmax.z, float(pmax.z));

	mmfini::pc ini;
	uint32 floatnbits = sizeof(double)*8;

	double3 const eigen_values = {eigen.GetEigenvalue(0), eigen.GetEigenvalue(1), eigen.GetEigenvalue(2)};
	double3 const
		eigen_vec0 = convert(eigen.GetEigenvector(0)),
		eigen_vec1 = convert(eigen.GetEigenvector(1)),
		eigen_vec2 = convert(eigen.GetEigenvector(2));

#define ASSIGN(field) bf::at_key<mmfini::field>(ini) = field
	ASSIGN(npoints);
	ASSIGN(floatnbits);
	ASSIGN(pmin);
	ASSIGN(pmax);
	ASSIGN(pcom);
	ASSIGN(p0);
	ASSIGN(eigen_values);
	ASSIGN(eigen_vec0);
	ASSIGN(eigen_vec1);
	ASSIGN(eigen_vec2);
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

	clog << "\r" << l << "lines (" << npoints << " points) in " << time(NULL) - start << "sec \n";
	return 0;
}