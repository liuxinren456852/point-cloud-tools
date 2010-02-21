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

class Grid {
	double find_i(double x) const { return (x - pmin.x) / res.x; }
	double find_j(double y) const { return (y - pmin.y) / res.y; }
public:
	uint64 nx, ny, ncells, npoints, ncellpts_max, ndatacells;
	double3 pmin, pmax, pcom, p0, res;

	#define LOAD(x) x = bf::at_key<mmfini::x>(ini)
	Grid(const mmfini::pc &ini, double xres0, double yres0) { 
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
		nx = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		ny = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		ncells = nx*ny;
	}
	Grid(const mmfini::pcsorted &ini) {
		LOAD(npoints);
		LOAD(pmin);
		LOAD(pmax);
		LOAD(pcom);
		LOAD(p0);
		LOAD(res);
		LOAD(ndatacells);
		LOAD(ncells);
		LOAD(ncellpts_max);

		nx = (uint64)((pmax.x-pmin.x)/res.x) + 1;
		ny = (uint64)((pmax.y-pmin.y)/res.y) + 1;
		const auto& dim = bf::at_key<mmfini::dim>(ini);
		BOOST_ASSERT(nx == dim.x && ny == dim.y);
		BOOST_ASSERT(ncells == nx*ny);
	}
	#undef LOAD

	uint64 index_encode (uint64 i, uint64 j) const { return i * ny + j;}
	void index_decode (uint64 ij, uint64 &i, uint64 &j) const { i = ij/ny; j = ij - i*ny;}
	uint64 index_encode(double x, double y) const { return index_encode(uint64(find_i(x)), uint64(find_j(y)));}

//idx = find(S==1);
//[ii,jj,kk] = ind2sub(size(S),idx) 
};

