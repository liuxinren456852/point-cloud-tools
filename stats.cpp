#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"
extern void calculate_cell_stats(const float3 points[], const uint64 n, Stats &stats) {
	double stdev = 0., skew = 0., kurt = 0.;
	stats.pzmax.z = -FLT_MAX;
	stats.pzmin.z = FLT_MAX;
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
	wstring filename(filename_ext), ext;
	const wstring::size_type beginext = filename_ext.rfind('.');
	if(beginext!=wstring::npos) {
		filename.assign(filename_ext, 0, beginext);
		ext.assign(filename_ext, beginext, wstring::npos);
	}

	// Open description file
	mmfini::pcsorted ini;
	if(!mmfini::load(ini, mmfini::filename(filename).c_str())) {
		cerr << "Cannot read the description file\n";
		return 1;
	}
	// Initialize the grid
	Grid grid(ini);
	// Validate parameters
	if(grid.floatnbits != 32) {
		clog << "The mmf appears to be unsorted (floatnbits != 32)\n";
		return 1;
	}
	// Open mmf files
	MMF<float3> points(filename_ext);
	MMF<Cell> cells(filename + _T("_cells") + ext);
	MMF<Stats> stats(filename + _T("_stats") + ext);
	if(	stats.open(bi::read_write)
		|| points.open()
		|| cells.open()
		)
		return 1;
	// Validate sizes
	BOOST_ASSERT(cells.size() == grid.ndatacells);
	BOOST_ASSERT(cells.size() == stats.size());
	clog << "\ncollecting statistics in data cells\n";
	clog << "assuming valid pcom data in the _stats file\n";
	clog << "using ";
	if(NUM_THREADS)
		omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
{
	#pragma omp master
	{
		clog << omp_get_num_threads() << " threads\n";
	}
	#pragma omp for schedule(dynamic, 1)
	for(auto ij = 0; ij < cells.size(); ij++) {
		const auto& cell = cells[ij];
		const float3 *const points_inthiscell = points.begin() + cell.startpt;
		calculate_cell_stats(points_inthiscell, cell.npoints, stats[ij]);
//		if(ij%1000 == 0)
//			clog << "\n" << ij << "    " << cell.startpt << "                ";
	}
} // #pragma omp parallel
	clog << "\r" << cells.size() << " cells\nstats() completed in " << time(NULL) - start << "sec\n";
	return err;
}