#include "stdafx.h"
#include "define.h"
#include "pctools.h"
#include "mmf.h"
#include "mmfini.h"
#include "grid.h"

static inline ostream& operator <<(ostream &os, const float3 &p) {
	os << p.x << ", " << p.y << ", " << p.z;
	return os;
}

int export_stats(const _TCHAR* pcmmf_filename, const uint64 min_npoints_per_cell)
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
	MMF<Cell> cells;
	MMF<Stats> stats;
	MMF<Stats> stats_detr;
	if(	stats.open(filename + _T("_stats") + ext)
		|| cells.open(filename + _T("_cells") + ext)
		|| stats_detr.open(filename + _T("_stats_detr") + ext)
		)
		return 1;
	// Validate sizes
	BOOST_ASSERT(grid.ndatacells == cells.size());
	BOOST_ASSERT(grid.ndatacells == stats.size());
	BOOST_ASSERT(grid.ndatacells == stats_detr.size());
	clog << "\nwriting the grid statistics... ";

	ofstream
		outmax((filename +_T("_zmax.txt")).c_str()),
		outmin((filename +_T("_zmin.txt")).c_str()),
		outcentroid((filename +_T("_centroid.txt")).c_str()),
		outstat((filename +_T("_stats.txt")).c_str()),
		outstat_underpopulated((filename+_T("_stats_underpopulated.txt")).c_str());

	// Prepare for output printing
	// Set Precision to print properly
	// Print the header
	//	outmax.precision(3);
	//	outmin.precision(3);
	//	outstat.precision(3);
	outmax.flags(ios::right + ios::fixed);
	outcentroid.flags(ios::right + ios::fixed);
	outmin.flags(ios::right + ios::fixed);
	outstat.flags(ios::right + ios::fixed);
	outstat_underpopulated.flags(ios::right + ios::fixed);
	outstat << "x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_detr,ku,ku_detr,zmean_d,n"<< endl;
	outstat_underpopulated << "x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_detr,ku,ku_detr,zmean_d,n"<< endl;
	outmin << "x,y,zmin" << endl;
	outmax << "x,y,zmax" << endl;
	outcentroid << "x,y,zmean" << endl;

	for(auto ij = 0; ij < cells.size(); ij++) {
		const auto& cell = cells[ij];
		const auto& s = stats[ij];
		const auto& sdetr = stats_detr[ij];
		uint64 i, j;
		grid.index_decode(cell.gridij, i, j);


		// Bring it back to the original system of coordinates
		const double xloc = grid.pmin.x+grid.res.x*i+grid.res.x/2   + grid.p0.x;
		const double yloc = grid.pmin.y+grid.res.y*j+grid.res.y/2   + grid.p0.y;
		const double3 pcom = s.pcom + grid.p0;
		const double3 pzmin = s.pzmin + grid.p0;
		const double3 pzmax = s.pzmax + grid.p0;

		uint64 nunderpopulatedcells = 0;

		// normally = zmean elevation at the center of the cell
		auto write = [&xloc, &yloc, &pcom, &pzmin, &pzmax, &s, &sdetr, &cell](ofstream &out) {
			out<<xloc<<","<<yloc<<","<<pcom.z << "," << pzmax.z<<","<<pzmin.z<<","<<pzmax.z-pzmin.z<<","<<s.stdev<<","<< sdetr.stdev<< ","<<s.skew <<","<<sdetr.skew<< "," <<s.kurt << "," <<sdetr.kurt<< "," << sdetr.pcom.z <<","<< cell.npoints<<endl;
		};

		if(cell.npoints < min_npoints_per_cell) {
			nunderpopulatedcells++;
			write(outstat_underpopulated);
//			outstat_underpopulated<<xloc<<","<<yloc<<","<<pcom.z << "," << pzmax.z<<","<<pzmin.z<<","<<pzmax.z-pzmin.z<<","<<s.stdev<<","<< sdetr.stdev<< ","<<s.skew <<","<<sdetr.skew<< "," <<s.kurt << "," <<sdetr.kurt<< "," << sdetr.pcom.z <<","<< cell.npoints<<endl;
			//			cell.stdev = cell.stdev_detrended = cell.skew = cell.kurt = NaN;
		} else {
			write(outstat);
			//	outstat << "x " << "y " << "zmean " << "zmax " << "zmin " << "range " << "stdev " << "stdev_detrended " << "zmean_detrended" << " n"<< endl;
			outmax<<pzmax << endl;
			outmin<<pzmin << endl;
			outcentroid<<pcom  << endl;
		}
	}
	clog << "done in " << time(NULL) - start << "sec\n\n";
	return err;
}