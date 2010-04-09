#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"

#include "hilbert.h"

inline void update_limits(const double3 &p, double3 &pmin, double3 &pmax) {
	if(p.x > pmax.x) pmax.x = p.x;
	if(p.y > pmax.y) pmax.y = p.y;
	if(p.z > pmax.z) pmax.z = p.z;

	if(p.x < pmin.x) pmin.x = p.x;
	if(p.y < pmin.y) pmin.y = p.y;
	if(p.z < pmin.z) pmin.z = p.z;
}

inline void update_limits(const Wm4::Vector3d &p, double3 &pmin, double3 &pmax) {
	if(p.X() > pmax.x) pmax.x = p.X();
	if(p.Y() > pmax.y) pmax.y = p.Y();
	if(p.Z() > pmax.z) pmax.z = p.Z();

	if(p.X() < pmin.x) pmin.x = p.X();
	if(p.Y() < pmin.y) pmin.y = p.Y();
	if(p.Z() < pmin.z) pmin.z = p.Z();
}
int sortintoxygrid(const _TCHAR* pcmmf_filename, uint32 xres_nom, uint32 yres_nom, uint32 denom, double normal0[3])
{
	int err = 0;
	clog << "\nSorting into a 2d grid";

	// File names
	const wstring filename_ext(pcmmf_filename);
	wstring filename(filename_ext), ext;
	const wstring::size_type beginext = filename_ext.rfind('.');
	if(beginext!=wstring::npos) {
		filename.assign(filename_ext, 0, beginext);
		ext.assign(filename_ext, beginext, wstring::npos);
	}

	mmfini::pc ini;
	if(!mmfini::load(ini, mmfini::filename(filename).c_str())) {
		cerr << "Cannot read the description file\n";
		return 1;
	}
	// Initialize grid
	Grid grid(ini, xres_nom, yres_nom, denom);
	clog << " with the resolution = " << grid.res << endl;

	if(grid.floatnbits != 64) {
		clog << "The mmf file appears sorted already. Transferring to changegridresolution()";
		return changegridresolution(pcmmf_filename, xres_nom, yres_nom, denom);
	}

	clog << "translating to -pcom\n";
	// Global detrending triggered by normal0 finished
	Wm4::Vector3d normal(0,0,1);
	if(normal0) {
		normal = normal0;
		clog << "rotating normal=" << normal << " --> [0,0,1]\n";
		normal.Normalize();
		Wm4::Vector3d const ez(0,0,1);
		Wm4::Quaterniond q;
		q.Align(normal, ez);
		Wm4::Vector3d slope(-ez + normal*normal.Z());
		slope.Normalize();
		clog << "rotating the 'slope', " << slope << ", --> [0,-1,0]\n";
		Wm4::Vector3d const neg_ey(0, -1, 0);
		Wm4::Vector3d const slope_after_1st_rotation = q.Rotate(slope);//-ez + normal*normal.Z());
		Wm4::Quaterniond q2;
		q2.Align(slope_after_1st_rotation, neg_ey);
		q = q2*q;
		clog << "rotation quaternion = " << q << endl;

		const auto pcom = grid.pcom;
		MMF<double3> points(filename_ext);
		if(points.open(bi::read_write))
			return 1;

		double3 pmin = {DBL_MAX,DBL_MAX,DBL_MAX}, pmax={-DBL_MAX,-DBL_MAX,-DBL_MAX};
		if(NUM_THREADS)
			omp_set_num_threads(NUM_THREADS);
		auto nteamthreads = omp_get_num_threads();

		#pragma omp parallel //reduction(+ : xx, xy, xz, yy, yz, zz)
		{
			double3 pmin_th = {DBL_MAX,DBL_MAX,DBL_MAX}, pmax_th={-DBL_MAX,-DBL_MAX,-DBL_MAX};
			#pragma omp for schedule(static, 1) nowait
			for(auto k = 0; k < grid.npoints; k++) {
				auto const p = points[k]-pcom;
				const Wm4::Vector3d pnew = q.Rotate(Wm4::Vector3d(as_array(p)));
				assign(points[k], pnew);
				update_limits(pnew, pmin_th, pmax_th);

			}
			#pragma omp critical 
			{
				update_limits(pmin_th, pmin, pmax);
				update_limits(pmax_th, pmin, pmax);
			}
			#pragma omp master
			nteamthreads = omp_get_num_threads();
		} // #pragma omp parallel finished

		clog << grid.npoints << " points, " << nteamthreads << " CPU threads\n";
		grid.translate(grid.pcom);
		grid.change_limits(pmin, pmax);

	} // global detrending triggered by normal0 finished


/*
	const uint64
		int_limit = std::numeric_limits<uint64>::max() - 1,
		ni_limit = uint64(  denom*(grid.pmax.x-grid.pmin.x)/xres_nom  ),
		nj_limit = uint64(  denom*(grid.pmax.y-grid.pmin.y)/yres_nom  ),
		npoints64 = bf::at_key<mmfini::npoints>(ini);
	clog << "is " << 8*sizeof(size_t) << "-bit address space enough? ...";
	if(npoints64 > int_limit) {
		clog << "no, too many points\n";
		return 1;
	}
	clog << "yes\n";
	if(ni_limit > int_limit || nj_limit > int_limit) {
		clog << "no, too many cells\n";
		return 1;
	}
	
	const auto npoints = uint64(npoints64);
*/

	clog << "number of cells at this resolution: "<< grid.size() << endl;
	
	// Memory-mapped files
	MMF<uint64>
		point_linked_list(filename + _T("_pll") + ext),
		head_point(filename + _T("_hp") + ext),
		ncellpts(filename + _T("_ncellpts") + ext),
		cellindex1(filename + _T("_cellindex1") + ext);
	MMF<double3> points(filename_ext);
	if(points.rename_and_open())
		return 1;
	MMF<float3>	points_sorted(filename_ext);
	if(point_linked_list.open(bi::read_write, grid.npoints+1)
		|| head_point.open(bi::read_write, grid.size())
		|| ncellpts.open(bi::read_write, grid.size())
		|| points_sorted.open(bi::read_write, grid.npoints)
		|| cellindex1.open(bi::read_write, grid.size())
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
	for(auto k = 0; k < grid.npoints; k++) {
//		const double3 p = points[k] - p0;
		const double3 p = points[k];
		auto const	ij = grid.find_cell(p);
		//ATOMIC ij ? impossible in OpenMP
		point_linked_list[k+1] = head_point[ij];
		head_point[ij] = k+1;
		ncellpts[ij]++;

		if(k%(1000000)==0) clog << "\r" << k;
	}
	clog << "\r" << grid.npoints << " points";
	clog << " in " << time(NULL) - start << "sec\n";











	// ndatacells
	uint64 ndatacells = 0;
	for(auto ij = 0; ij < ncellpts.size(); ij++) {
		if(ncellpts[ij] > 0)
			ndatacells++;
	}
	grid.ndatacells = ndatacells;












	// Now we know the number of non-empty data cells
	MMF<Cell> cells(filename + _T("_cells") + ext);
	if(	cells.open(bi::read_write, grid.ndatacells)	)
		return 1;


	// Data cells will hold pointers into the point cloud
	grid.ncellpts_max = 0;
	uint64 datacellindex = 0;
	auto pcell = &cells[0];
	uint64 startpt = 0;
//	for(auto ij = 0; ij < ncellpts.size(); ij++) {
#ifdef HILBERT
	clog << "Cells sorting: Hilbert order\n";
	const Hilbert2DSpaceFillingCurve sfc( grid.nx(),grid.ny() );
	for(uint64 h = 0; h < sfc.length(); h++) {
		uint64 i, j;
		sfc.decode(h, i, j);
		if(i < grid.nx() && j < grid.ny()) {
			const auto ij = grid.index_encode(i, j);

#else
	clog << "Cells sorting: XRowMajor\n";
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
	MMF<Stats> stats(filename + _T("_stats") + ext);
	if(	stats.open(bi::read_write, grid.ndatacells)	)
		return 1;
//	fill(stats.begin(), stats.end(), g_init_stats);



	start = time(0);
	clog << "writing the sorted point cloud... " << endl;
	clog << "using ";
	if(NUM_THREADS)
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
{
	// Gather the points into contiguous blocks and write the sorted binary file
	// PARALLEL
	vector<double3> points_inthiscell;
	points_inthiscell.reserve(grid.ncellpts_max);

	#pragma omp for schedule(dynamic, 1) nowait
	for(auto ij = 0; ij < ndatacells; ij++) {
		const auto& cell = cells[ij];
		const auto n = cell.npoints;
		points_inthiscell.resize(n);
		auto points_inthiscell_it = points_inthiscell.begin();
		// Stats
		auto& cellstats = stats[ij];
		//cellstats.pcom = {0,0,0}; // Assuming done
		float3 *points_inthiscell_save = points_sorted.begin() + cell.startpt;
		double3 pcom_inthiscell = {0,0,0};
		auto k1 = head_point[cell.gridij];
		//uint64 ncellpoints = 0;
		while(k1) {
			// Shift the point relative to grid.pcom
			auto const p = points[k1-1] - grid.pcom;
			pcom_inthiscell += p;
			*points_inthiscell_it++ = p;
//			if(vertical_project)
//				p.z = dot(p, normal/nz);

			// Next point
			k1 = point_linked_list[k1];
		}
		//BOOST_ASSERT(ncellpoints == cell.npoints);

		// Finalize pcom, note that it is relative to the grid.pcom
		pcom_inthiscell /= n;
		cellstats.pcom = convert(pcom_inthiscell);
		
		// We'll do some stats
		double stdev = 0., skew = 0., kurt = 0.;
		cellstats.pzmax.z = -FLT_MAX;
		cellstats.pzmin.z = FLT_MAX;

		for(auto k = 0; k < n; k++) {
			const auto p = points_inthiscell[k] - pcom_inthiscell;
			*points_inthiscell_save++ = convert(p);

			const double dz = p.z;// - cellstats.pcom.z;
			const double dz2 = dz * dz;
			stdev += dz2;
			skew += dz2 * dz;
			kurt += dz2 * dz2;

			if (p.z > cellstats.pzmax.z)
				cellstats.pzmax = convert(p);
			if (p.z < cellstats.pzmin.z) 
				cellstats.pzmin = convert(p);
		}
		// Finalize stats for this cell
		if(n > 1) {
			double const stdev2 = stdev / (n-1);
			double const stdev = sqrt(stdev2);
			cellstats.stdev = stdev;
			if(n > 3) {
				cellstats.skew = (double(n) / ((n-1)*(n-2))) * skew / (stdev2 * stdev);
				kurt = (double(n) * (n+1) / (n-1)) * kurt / (stdev2 * stdev2);
				kurt -= 3. * (n-1) * (n-1);
				kurt /= ((n-2)*(n-3));
				cellstats.kurt = kurt;
			}
		}

		points_inthiscell.clear();
	}
#pragma omp master
	{
		clog << omp_get_num_threads() << " threads\n";
	}
} // #pragma omp parallel

		// Let's not forget to update the grid as we translated each point and each cell's center
		grid.translate(grid.pcom);

		// Delete old unsorted point cloud and auxiliary files
		points.delete_file();
		point_linked_list.delete_file();
		head_point.delete_file();

		mmfini::pcsorted ini_sorted;
		#define ASSIGN(x) bf::at_key<mmfini::x>(ini_sorted) = grid.x
		bf::at_key<mmfini::npoints>(ini_sorted) = grid.npoints;
		bf::at_key<mmfini::floatnbits>(ini_sorted) = sizeof(float)*8;
		ASSIGN(pmin);
		ASSIGN(pmax);
		ASSIGN(pcom);
		ASSIGN(p0);
		bf::at_key<mmfini::eigen_values>(ini_sorted) = bf::at_key<mmfini::eigen_values>(ini);
		bf::at_key<mmfini::eigen_vec0>(ini_sorted) = bf::at_key<mmfini::eigen_vec0>(ini);
		bf::at_key<mmfini::eigen_vec1>(ini_sorted) = bf::at_key<mmfini::eigen_vec1>(ini);
		bf::at_key<mmfini::eigen_vec2>(ini_sorted) = bf::at_key<mmfini::eigen_vec2>(ini);
//		ASSIGN(eigen_values);
//		ASSIGN(eigen_vec0);
//		ASSIGN(eigen_vec1);
//		ASSIGN(eigen_vec2);
		double3 normal_save = *(double3*)(double*)normal;
		bf::at_key<mmfini::grid_plane_normal>(ini_sorted) = normal_save;
		//ASSIGN(grid_plane_normal);
		ASSIGN(res);
		ASSIGN(dim);
		//const Point3<uint64> dim = {grid.nx, grid.ny, 0};
		//bf::at_key<mmfini::dim>(ini_sorted) = dim;

		bf::at_key<mmfini::ncells>(ini_sorted) = grid.size();
//		ASSIGN(ncells);
		ASSIGN(ndatacells);
		ASSIGN(ncellpts_max);
		#undef ASSIGN
		if(!mmfini::save(ini_sorted, (mmfini::filename(filename)).c_str())) {
			cerr << "Cannot save .ini file";
			return 1;
		}
	clog << "sortintoxygrid() completed in " << time(NULL) - start << "sec\n";
	return err;
}


/*
Wm4::Vector3d const
normal = eigen.GetEigenvector(0),
ez(0,0,1),
ey(0,1,0);

Wm4::Quaterniond q1,q2;
q1.Align(normal, ez);

Wm4::Vector3d slope = q1.Rotate(-ez + normal*normal.Z());
slope.Normalize();
q2.Align(slope, -ey);

Wm4::Quaterniond const q = q2*q1;
//	clog << "q = " << q.W() << ' ' << q.X() << ' ' << q.Y() << ' ' << q.Z() << endl;
clog << "q = " << q << endl;

Wm4::Vector3d const
normal_new = q.Rotate(normal),
slope_new = q.Rotate(-ez + normal*normal.Z());
*/