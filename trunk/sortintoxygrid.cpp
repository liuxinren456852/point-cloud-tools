#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"

#define HILBERT
#include "hilbert.h"

int sortintoxygrid(const _TCHAR* file, double xres, double yres)
{
	int err = 0;

	const wstring filename_ext(file);
	const wstring::size_type beginext = filename_ext.rfind('.');
	const wstring filename(filename_ext, 0, beginext);
	const wstring ext(filename_ext, beginext);

	mmfini::pc ini;
	if(!mmfini::load(ini, mmfini::filename(filename).c_str())) {
		cerr << "Cannot read the description file\n";
		return 1;
	}
	// Rename the current unsorted file 
	// The .ini filename stays the same and will be overwritten
	const wstring unsortedfile = wstring(file) + _T("_unsorted");
	if(	_wrename(file, unsortedfile.c_str())) {
		cerr << "Cannot rename the file\n";
		return 1;
	}

	// Initialize grid
	Grid grid(ini, xres, yres);
	
	const uint64
		int_limit = std::numeric_limits<uint64>::max() - 1,
		ni_limit = uint64((grid.pmax.x-grid.pmin.x)/xres),
		nj_limit = uint64((grid.pmax.y-grid.pmin.y)/yres),
		npoints64 = bf::at_key<mmfini::npoints>(ini);
	clog << "is " << 8*sizeof(size_t) << "-bit address space enough? ...";
	if(npoints64 > int_limit) {
		clog << "no, too many points\n";
		return 1;
	}
	if(ni_limit > int_limit || nj_limit > int_limit) {
		clog << "no, too many cells\n";
		return 1;
	}
	clog << "yes\n";
 	const auto npoints = uint64(npoints64);

	clog << "number of cells at this resolution: "<< grid.size << endl;
	
	// Memory-mapped files
	MMF<uint64> point_linked_list, head_point, ncellpts;
	MMF<double3> points;
	MMF<float3>	points_sorted;
	MMF<uint64>	cellindex1;
	if(point_linked_list.open(filename + _T("_pll") + ext, bi::read_write, npoints+1)
		|| head_point.open(filename + _T("_hp") + ext, bi::read_write, grid.memsize)
		|| ncellpts.open(filename + _T("_ncellpts") + ext, bi::read_write, grid.memsize)
		|| points.open(unsortedfile)
		|| points_sorted.open(filename_ext, bi::read_write, npoints)
		|| cellindex1.open(filename + _T("_cellindex1") + ext, bi::read_write, grid.memsize)
		)
		return 1;


	// "The _chsize_s function extends or truncates the file associated with fd to the length specified by size.
	// "...Null characters ('\0') are appended if the file is extended."
	// So it should be filled with 0 already
	//fill(grid.begin(),grid.end(),g_init_cell);
	//	fill(stats, stats+grid.size(),g_init_stats);






















	// 
	// Read a point from unsorted binary point cloud
	// Find the corresponding grid cell
	// Update the linked list
	// And we can calculate already some cell statistics: zmean, pmin, and pmax
	// PARALLEL? No. Because of the nature of linked list building
	clog << "indexing points into the cells... \n";
	time_t start = time(NULL);
	const double3 p0 = grid.p0;
	for(auto k = 0; k < npoints; k++) {
		const double3 p = points[k] - p0;
		auto const	ij = grid.find_cell(p.x, p.y);
		//ATOMIC ij ? impossible in OpenMP
		point_linked_list[k+1] = head_point[ij];
		head_point[ij] = k+1;
		ncellpts[ij]++;

		if(k%(1000000)==0) clog << "\r" << k;
	}
	clog << "\r" << npoints << " points";
	clog << " in " << time(NULL) - start << "sec\n";











	// ndatacells
	uint64 ndatacells = 0;
	for(auto ij = 0; ij < ncellpts.size(); ij++) {
		if(ncellpts[ij] > 0)
			ndatacells++;
	}
	grid.ndatacells = ndatacells;












	// Now we know the number of non-empty data cells
	MMF<Cell> cells;
	if(	cells.open(filename + _T("_cells") + ext, bi::read_write, grid.ndatacells)	)
		return 1;


	// Data cells will hold pointers into the point cloud
	grid.ncellpts_max = 0;
	uint64 datacellindex = 0;
	auto pcell = &cells[0];
	uint64 startpt = 0;
//	for(auto ij = 0; ij < ncellpts.size(); ij++) {
#ifdef HILBERT
	Hilbert2DSpaceFillingCurve sfc( grid.nx(),grid.ny() );
	for(uint64 h = 0; h < sfc.length(); h++) {
		uint64 i, j;
		sfc.decode(h, i, j);
		if(i < grid.nx() && j < grid.ny()) {
			const auto ij = grid.index_encode(i, j);

#else
	{

		for(auto ij = 0; ij < ncellpts.size(); ij++) {
#endif
			auto const npts = ncellpts[ij];
			if(npts > 0) {
				cellindex1[ij] = 1 + datacellindex++;
				auto& cell = *pcell++; 
				cell.gridij = ij;
				cell.startpt = startpt;
				cell.npoints = npts;
				startpt+=npts;
				if(npts > grid.ncellpts_max) grid.ncellpts_max = npts;
			}
		}
	}
	ncellpts.delete_file();
	BOOST_ASSERT(datacellindex == ndatacells);

	clog << grid.ndatacells << " data cells\n";
	clog << grid.ncellpts_max << " maximum number of points per cell\n";




	// Sorting the point cloud

	// We will also start calculating cell statistics
	MMF<Stats> stats;
	if(	stats.open(filename + _T("_stats") + ext, bi::read_write, grid.ndatacells)	)
		return 1;
	fill(stats.begin(), stats.end(), g_init_stats);



	start = time(0);
	clog << "writing the sorted point cloud... " << endl;
	//	omp_set_num_threads(1);
	//omp_set_num_threads(omp_get_num_procs());
	// Let's try 10 threads, probably more than 
	clog << "using ";
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
{
	#pragma omp master
	{
		clog << omp_get_num_threads() << " threads\n";
	}
	// Gather the points into contiguous blocks and write the sorted binary file
	// PARALLEL
	#pragma omp for
	for(auto ij = 0; ij < ndatacells; ij++) {
		const auto& cell = cells[ij];
		// Stats
		auto& cellstats = stats[ij];
		//cellstats.pcom = {0,0,0}; // Assuming done
		float3 *points_inthiscell = points_sorted.begin() + cell.startpt;
		double3 pcom = {0,0,0};
		auto k1 = head_point[cell.gridij];
		while(k1) {
			// Shift the point to grid.p0 (old pcom)
			auto const p = points[k1-1] - grid.p0;
			// Copy the point truncated to float
			*points_inthiscell++ = convert(p);
			// Some cell stats
			pcom += p;
			// Next point
			k1 = point_linked_list[k1];
		}
		// Finalize stats
		pcom /= cell.npoints;
		cellstats.pcom = convert(pcom);
	}
} // #pragma omp parallel

		// Delete old unsorted point cloud and auxiliary files
		points.delete_file();
		point_linked_list.delete_file();
		head_point.delete_file();

		clog << "\r" << grid.ndatacells << " cells";
		clog << " in " << time(NULL) - start << "sec\n";

		mmfini::pcsorted ini_sorted;
		#define ASSIGN(x) bf::at_key<mmfini::x>(ini_sorted) = grid.x
		bf::at_key<mmfini::npoints>(ini_sorted) = npoints;
		bf::at_key<mmfini::floatnbits>(ini_sorted) = sizeof(float)*8;
		ASSIGN(pmin);
		ASSIGN(pmax);
		ASSIGN(pcom);
		ASSIGN(p0);
		ASSIGN(res);
		ASSIGN(dim);
		//const Point3<uint64> dim = {grid.nx, grid.ny, 0};
		//bf::at_key<mmfini::dim>(ini_sorted) = dim;

		bf::at_key<mmfini::ncells>(ini_sorted) = grid.size;
//		ASSIGN(ncells);
		ASSIGN(ndatacells);
		ASSIGN(ncellpts_max);
		#undef ASSIGN
		if(!mmfini::save(ini_sorted, (mmfini::filename(filename)).c_str())) {
			cerr << "Cannot save .ini file";
			err = 1;
		}

	clog << "sortintoxygrid() completed with code " << err << endl << endl;
	return err;
}