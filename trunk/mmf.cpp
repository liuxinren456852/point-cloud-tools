#include "stdafx.h"
#include "define.h"
namespace bi = boost::interprocess;
int init_mmf(bi::file_mapping &fm, bi::mapped_region &mr, const _TCHAR* file, bi::mode_t mode, int64 resize_bytes) {
	const wstring wfilename(file);
	if(resize_bytes > 0) // Resize the file
	{
		FILE *file = _tfopen(wfilename.c_str(), _T("wbS"));
		if(file == NULL) {
			cerr << "Cannot open a mmf file for writing " << endl; return 1;}
		if(_chsize_s(_fileno(file), resize_bytes) != 0) {//EBADF 
			cerr << "Can't resize the output file to" << resize_bytes << endl; return 1; }
		fclose(file);
	}
	// Current limitation of boost::interprocess, no wchar strings
	const string filename(wfilename.begin(), wfilename.end());
	try {
		bi::file_mapping fm1(filename.c_str(),	mode);
		fm.swap(fm1);
		bi::mapped_region mr1(fm,				mode);
		mr.swap(mr1);
	}
	catch(bi::interprocess_exception &ex){
		cout << ex.what() << endl;
		return 1;
	}
	if(resize_bytes > 0 && mr.get_size() != resize_bytes) {
		cerr << "File size is wrong:" << mr.get_size() << " against " << resize_bytes << endl; return 1;
	}
	return 0;
}