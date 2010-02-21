#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"
inline void calculate_cell_stats(const float3 points[], const uint64 n, Stats &stats) {
	double stdev = 0., skew = 0., kurt = 0.;
	for(auto k = 0; k < n; k++) {
		auto p = points[k];
		const double dz = p.z - stats.pcom.z;
		const double dz2 = dz * dz;
		stdev += dz2;
		skew += dz2 * dz;
		kurt += dz2 * dz2;

		if (p.z > stats.pzmax.z)
			stats.pzmax = p;
		if (p.z < stats.pzmin.z) 
			stats.pzmin = p;
	}
	// Finalize stats for this cell
	if(n > 1) {
		double const stdev2 = stdev / (n-1);
		double const stdev = sqrt(stdev2);
		stats.stdev = stdev;
		if(n > 3) {
			stats.skew = (double(n) / ((n-1)*(n-2))) * skew / (stdev2 * stdev);
			kurt = (double(n) * (n+1) / (n-1)) * kurt / (stdev2 * stdev2);
			kurt -= 3. * (n-1) * (n-1);
			kurt /= ((n-2)*(n-3));
			stats.kurt = kurt;
		}
	}
}
int stats(const _TCHAR* pcmmf_filename)
{
	int err = 0;
	time_t start = time(NULL);

	// File names
	const wstring filename_ext(pcmmf_filename);
	const wstring::size_type beginext = filename_ext.rfind('.');
	const wstring filename(filename_ext, 0, beginext);
	const wstring ext(filename_ext, beginext);
	// Open description file
	mmfini::pcsorted ini;
	if(!mmfini::load(ini, mmfini::filename(filename).c_str())) {
		cerr << "Cannot read the description file\n";
		return 1;
	}
	// Initialize the grid
	Grid grid(ini);
	// Validate parameters
	BOOST_ASSERT(bf::at_key<mmfini::floatnbits>(ini) == sizeof(float)*8);
	// Open mmf files
	MMF<float3> points;
	MMF<Cell> cells;
	MMF<Stats> stats;
	if(	stats.open(filename + _T("_stats") + ext, bi::read_write)
		|| points.open(filename_ext)
		|| cells.open(filename + _T("_cells") + ext)
		)
		return 1;
	// Validate sizes
	BOOST_ASSERT(cells.size() == grid.ndatacells);
	BOOST_ASSERT(cells.size() == stats.size());
	clog << "\ncollecting statistics in data cells\n";
	clog << "assuming valid pcom data in the _stats file\n";
	clog << "using ";
	omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
{
	#pragma omp master
	{
		clog << omp_get_num_threads() << " threads\n";
	}
	#pragma omp for
	for(auto ij = 0; ij < cells.size(); ij++) {
		const auto& cell = cells[ij];
		const float3 *const points_inthiscell = points.begin() + cell.startpt;
		calculate_cell_stats(points_inthiscell, cell.npoints, stats[ij]);
	}
} // #pragma omp parallel
	clog << "stats() completed in " << time(NULL) - start << "sec\n";
	return err;
}