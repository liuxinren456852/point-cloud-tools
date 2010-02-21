#pragma once
namespace bf = boost::fusion;
namespace mmfini
{
	#define FIELD( x, t ) struct x {typedef t type; typedef bf::pair<x, t> pair; static const char* name() {return #x;};}
	FIELD(npoints, uint64);
	FIELD(floatnbits, uint32);
// 	FIELD(xmin, float64);
// 	FIELD(xmax, float64);
// 	FIELD(ymin, float64);
// 	FIELD(ymax, float64);
// 	FIELD(zmin, float64);
// 	FIELD(zmax, float64);
	FIELD(ncells, uint64);
	FIELD(ndatacells, uint64);
	FIELD(ncellpts_max, uint64);
	FIELD(res, double3);
// 	FIELD(yres, float64);
// 	FIELD(xmean, float64);
// 	FIELD(ymean, float64);
// 	FIELD(zmean, float64);
// 	FIELD(x0, float64);
// 	FIELD(y0, float64);
// 	FIELD(z0, float64);
	FIELD(pmin, double3);
	FIELD(pmax, double3);
	FIELD(pcom, double3);
	FIELD(p0, double3);
	FIELD(dim, Point3<uint64>);
//	FIELD(ny, uint64);
//	FIELD(nz, uint64);
	#undef FIELD

	typedef bf::map	<
		npoints::pair,
		floatnbits::pair,
// 		xmin::pair,
// 		xmax::pair,
// 		ymin::pair,
// 		ymax::pair,
// 		zmin::pair,
// 		zmax::pair,
// 		xmean::pair,
// 		ymean::pair,
// 		zmean::pair
		pmin::pair,
		pmax::pair,
		pcom::pair,
		p0::pair
	> pc;
	typedef bf::map	<
		npoints::pair,
		floatnbits::pair,
		pmin::pair,
		pmax::pair,
		pcom::pair,
		p0::pair,
		res::pair,
		dim::pair,
		ncells::pair,
		ndatacells::pair,
		ncellpts_max::pair	> pcsorted;
// 	typedef bf::map	<
// 		nx::pair,
// 		ny::pair,
// 		floatnbits::pair,
// 		x0::pair,
// 		y0::pair,
// 		z0::pair,
// 		xmin::pair,
// 		xmax::pair,
// 		ymin::pair,
// 		ymax::pair,
// 		zmin::pair,
// 		zmax::pair,
// 		xres::pair,
// 		yres::pair	> grid;
// 	typedef bf::map	<
// 		count::pair,
// 		floatnbits::pair,
// 		nx::pair,
// 		ny::pair,
// 		nz::pair,
// 		xmin::pair,
// 		xmax::pair	> grid3d;
	typedef bf::map	<
//		nx::pair,
		dim::pair	> grid;
// 	typedef bf::map	<
// 		nx::pair,
// 		ny::pair,
// 		nz::pair	> matrix3D;

	struct save_field {
		ofstream &file;
		template <typename Pair>
		void operator()(Pair const& data) const {
//			file << typeid(Pair::first_type).name()+fields::typename_offset << " = " << data.second << ";\n";
			file << Pair::first_type::name() << " = " << data.second << ";\n";
		}
		save_field(ofstream &file) : file(file) {};
	};
	struct load_field {
		ifstream &file;
		template <typename Pair>
		void operator()(Pair & data) const {
			file >> Pair::first_type::name() >> ws >> '=' >> data.second >> ws >>';';
			//clog << data.second <<endl;
		}
		load_field(ifstream &file) : file(file) {};
	};

	template <typename Stuff>
	bool save(Stuff const& stuff, const _TCHAR* name) {
		ofstream file(name);
		//file.flags(ios::fixed);
		file.precision(std::numeric_limits<double>::digits10);
		bf::for_each(stuff, save_field(file));
		return file.good();
	};

	template <typename Stuff>
	bool load(Stuff & stuff, const _TCHAR* name) {
		ifstream file(name);
		bf::for_each(stuff, load_field(file));
		return file.good();
	};

	static wstring filename(const _TCHAR* mmf) {return wstring(mmf) + _T(".ini");}
	static wstring filename(const wstring &mmf) {return mmf + _T(".ini");}
};
/*
// Example usage
	mmfini::xyz x, x1;
	bf::at_key<mmfini::floatnbits>(x) = 64;
	bf::at_key<mmfini::xmin>(x) = 4.1;
	bf::at_key<mmfini::xmax>(x) = 5.1;
	uint32 floatnbits = bf::at_key<mmfini::floatnbits>(x);
	mmfini::save(x, _T("test.txt"));
	if(mmfini::load(x1, _T("test.txt"))) {
		uint32 floatnbits = at_key<mmfini::floatnbits>(x1);
	}
	mmfini::save(x1, _T("test.txt"));
*/