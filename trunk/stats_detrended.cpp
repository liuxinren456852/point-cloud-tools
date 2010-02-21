#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"
void calculate_cell_stats(const float3 points[], const uint64 n, Stats &stats);

int stats_detrended(const _TCHAR* pcmmf_filename, const double times_zmin, const double times_zmean, const double times_stdev)
{
	// The ground level is taken to be at the center of the cell;
	// it is a linear combination of cell min, average elevations and the stdev
	// normally = zmean elevation at the center of the cell
	auto ground_level = [&times_zmin, &times_zmean, &times_stdev](const Stats &cell_stats) {
		return times_zmin*double(cell_stats.pzmin.z) + times_zmean*double(cell_stats.pcom.z) + times_stdev*double(cell_stats.stdev);
	};




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
	MMF<Stats> stats_detr;
	MMF<uint64>	cellindex1;
	if(	stats.open(filename + _T("_stats") + ext)
		|| points.open(filename_ext)
		|| cells.open(filename + _T("_cells") + ext)
		|| stats_detr.open(filename + _T("_stats_detr") + ext, bi::read_write, grid.ndatacells)
		|| cellindex1.open(filename + _T("_cellindex1") + ext)
		)
		return 1;
	// Validate sizes
	BOOST_ASSERT(grid.ndatacells == cells.size());
	BOOST_ASSERT(grid.ndatacells == stats.size());
	BOOST_ASSERT(grid.ndatacells == stats_detr.size());
	BOOST_ASSERT(grid.ncells == cellindex1.size());
	clog << "\ncollecting 'detrended' statistics in data cells\n";
	clog << "using ";
	omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
{
	#pragma omp master
	{
		clog << omp_get_num_threads() << " threads\n";
 	}

	// We will need a per-thread buffer for detrended points in each cell
	vector<float3> points_inthiscell_detrended;
	points_inthiscell_detrended.reserve(grid.ncellpts_max);

	#pragma omp for
	for(auto ij = 0; ij < cells.size(); ij++) {
		const auto& cell = cells[ij];
		const auto& cellstats = stats[ij];
		auto& cellstats_detr = stats_detr[ij];
		uint64 i, j;
		grid.index_decode(cell.gridij, i, j);
		const double xloc=grid.pmin.x+grid.res.x*i+grid.res.x/2;
		const double yloc=grid.pmin.y+grid.res.y*j+grid.res.y/2;
		auto const
			iplus1 = min<uint64>(i + 1, grid.nx-1),
			iminus1 = max<uint64>(i, 1) - 1,
			jplus1 = min<uint64>(j + 1, grid.ny-1),
			jminus1 = max<uint64>(j, 1) - 1;




		const float3 *const points_inthiscell = points.begin() + cell.startpt;
		double3 pcom = {0,0,0};
		const auto n = cell.npoints;
		points_inthiscell_detrended.resize(n);
		for(auto k = 0; k < n; k++) {
			auto p = points_inthiscell[k];


			auto const
				ix = p.x > xloc ? iplus1 : iminus1,
				jy = p.y > yloc ? jplus1 : jminus1;
			const auto neighbour_x_cell_index1 = cellindex1[grid.index_encode(ix, j)];
			const auto neighbour_y_cell_index1 = cellindex1[grid.index_encode(i, jy)];
			double const
				z0 = ground_level(cellstats),
				zx = neighbour_x_cell_index1 > 0 ? ground_level(stats[neighbour_x_cell_index1-1]) : z0,
				zy = neighbour_y_cell_index1 > 0 ? ground_level(stats[neighbour_y_cell_index1-1]) : z0;
			//Detrend
			p.z -= (z0 + (zx - z0) * abs(p.x - xloc)/grid.res.x + (zy - z0) * abs(p.y - yloc)/grid.res.y);

			points_inthiscell_detrended[k] = p;
			pcom += convert(p);
		}
		pcom /= n;
		cellstats_detr.pcom = convert(pcom);
		calculate_cell_stats(&points_inthiscell_detrended[0], n, cellstats_detr);
		points_inthiscell_detrended.clear();
	}
} // #pragma omp parallel
	clog << "stats_detrended() completed in " << time(NULL) - start << "sec\n";
	return err;
}