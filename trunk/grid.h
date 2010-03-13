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
Cell const g_init_cell = {0, 0, 0};
Stats const g_init_stats = {
	{0.F, 0.F, 0.F}, {0.F, 0.F, FLT_MAX}, {0.F, 0.F, -FLT_MAX},
	0.F, 0.F, 0.F
};

class XRowMajor {
protected:
	void init(const Point3<uint64> &dim) {
		BOOST_ASSERT(dim.z == 1); // no 3D
		memnx = dim.x;
		memny = dim.y;
		memsize = memnx*memny;
	}
public:
	uint64 memnx, memny; // the actual grids in memory, == dim of the bounding box
	uint64 memsize; // this is the actual number of cells in memory, == nx*ny
	uint64 index_encode (uint64 i, uint64 j) const { return i * memny + j;}
	const uint64 & nx () const { return memnx;}
	const uint64 & ny () const { return memny;}
	void index_decode (uint64 ij, uint64 &i, uint64 &j) const { i = ij/memny; j = ij - i*memny;}
};
template<typename MemoryLayout>
class GridT : public MemoryLayout {
	double find_i(double x) const { return (x - pmin.x) / res.x; }
	double find_j(double y) const { return (y - pmin.y) / res.y; }
public:
	uint64 size, npoints, ncellpts_max, ndatacells;
	Point3<uint64> dim; // index dimensions from the bounding box
	double3 pmin, pmax, pcom, p0, res;

	#define LOAD(x) x = bf::at_key<mmfini::x>(ini)
	GridT(const mmfini::pc &ini, double xres0, double yres0) { 
		LOAD(pmin);
		LOAD(pmax);
		LOAD(pcom);
		LOAD(p0);

		p0 += pcom;
		pmin -= pcom;
		pmax -= pcom;
		pcom -= pcom;

		res.x = xres0;
		res.y = yres0;
		res.z = 0;

		dim.x = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		dim.y = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		dim.z = 1;

		size = dim.x * dim.y * dim.z;
		MemoryLayout::init(dim);
		//nx = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		//ny = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		//ncells = nx*ny;
	}
	GridT(const mmfini::pcsorted &ini) {
		LOAD(npoints);
		LOAD(pmin);
		LOAD(pmax);
		LOAD(pcom);
		LOAD(p0);
		LOAD(res);
		LOAD(ndatacells);
		//LOAD(ncells);
		LOAD(ncellpts_max);

		dim.x = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		dim.y = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		dim.z = 1;
		size = dim.x * dim.y * dim.z;
		MemoryLayout::init(dim);
		const auto& dim0 = bf::at_key<mmfini::dim>(ini);
		//const auto& ncells = bf::at_key<mmfini::ncells>(ini);
		BOOST_ASSERT(dim.x == dim0.x && dim.y == dim0.y && dim.z == dim0.z);
		//BOOST_ASSERT(ncells == dim.x*dim.y*dim.z);
	}
	#undef LOAD

	uint64 find_cell(double x, double y) const { return MemoryLayout::index_encode(uint64(find_i(x)), uint64(find_j(y)));}

//idx = find(S==1);
//[ii,jj,kk] = ind2sub(size(S),idx) 
};

typedef GridT<XRowMajor> Grid;
