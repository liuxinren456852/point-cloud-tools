#pragma once


struct Cell {
	uint64 gridij;
	uint64 npoints;
	uint64 startpt;
};
struct Stats {
	float3 pcom, pzmin, pzmax;
	float stdev, skew, kurt;
};
// Cell const g_init_cell = {0, 0, 0};
// Stats const g_init_stats = {
// 	{0.F, 0.F, 0.F}, {0.F, 0.F, FLT_MAX}, {0.F, 0.F, -FLT_MAX},
// 	0.F, 0.F, 0.F
// };


class Grid {
public:
	const uint64 & nx () const { return dim.x;}
	const uint64 & ny () const { return dim.y;}
	const uint64 size () const { return dim.x*dim.y;}
	void index_decode (const uint64 ij, uint64 &i, uint64 &j) const { i = ij/ny(); j = ij - i*ny();}
	uint64 index_encode (uint64 i, uint64 j) const { return i * ny() + j;}

	uint64 npoints, ncellpts_max, ndatacells;
	Point3<uint64> dim; // index dimensions from the bounding box
	double3 pmin, pmax, pcom, p0;
	Tuple4<uint32> res;
	uint8 floatnbits;

	#define LOAD(x) x = bf::at_key<mmfini::x>(ini)
	Grid(const mmfini::pc &ini, uint32 xres_nom, uint32 yres_nom, uint32 denom) { 
		LOAD(npoints);
		LOAD(floatnbits);
		LOAD(pmin);
		LOAD(pmax);
		LOAD(pcom);
		LOAD(p0);
		ncellpts_max = 0;
		ndatacells = 0;

		//p0 += pcom;
		//pmin -= pcom;
		//pmax -= pcom;
		//pcom -= pcom;

		res.x = xres_nom;
		res.y = yres_nom;
		res.z = 0;
		res.w = denom;

		dim.x = find_i(pmax.x) + 1;
		dim.y = find_j(pmax.y) + 1;
		dim.z = 1;

		//MemoryLayout::init(dim);
		//nx = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		//ny = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		//ncells = nx*ny;
	}
	Grid(const mmfini::pcsorted &ini) {
		LOAD(npoints);
		LOAD(floatnbits);
		LOAD(pmin);
		LOAD(pmax);
		LOAD(pcom);
		LOAD(p0);
		LOAD(res);
		LOAD(ndatacells);
		//LOAD(ncells);
		LOAD(ncellpts_max);

		dim.x = find_i(pmax.x) + 1;
		dim.y = find_j(pmax.y) + 1;
		dim.z = 1;

		//MemoryLayout::init(dim);
		const auto& dim0 = bf::at_key<mmfini::dim>(ini);
		//const auto& ncells = bf::at_key<mmfini::ncells>(ini);
		BOOST_ASSERT(dim.x == dim0.x && dim.y == dim0.y && dim.z == dim0.z);
		//BOOST_ASSERT(ncells == dim.x*dim.y*dim.z);
	}
	#undef LOAD

	void translate(const double3 &p) {
		p0 += p;
		pmin -= p;
		pmax -= p;
		pcom -= p;
	}
	void change_limits(const double3 &pmin_new, const double3 &pmax_new) {
		pmin = pmin_new;
		pmax = pmax_new;
		dim.x = find_i(pmax.x) + 1;
		dim.y = find_j(pmax.y) + 1;
		dim.z = 1;
	}
	uint64 find_i(double x) const { return res.w*(x - pmin.x) / res.x; }
	uint64 find_j(double y) const { return res.w*(y - pmin.y) / res.y; }
//	uint64 find_cell(double x, double y) const { return MemoryLayout::index_encode(find_i(x), find_j(y));}
	template<typename T>
	uint64 find_cell(const Point3<T> &p) const { return index_encode(find_i(p.x), find_j(p.y));}

//idx = find(S==1);
//[ii,jj,kk] = ind2sub(size(S),idx) 
};

