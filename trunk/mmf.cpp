#include "stdafx.h"
#include "define.h"
namespace bi = boost::interprocess;

int resize(const wstring &filename, int64 resize_bytes)
{
	FILE *file = _tfopen(filename.c_str(), _T("abS"));
	if(file == NULL) {
		wcerr << "Cannot open " << filename <<"for writing\n"; return 1;}
	if(_chsize_s(_fileno(file), resize_bytes) != 0) {//EBADF 
		wcerr << "Cannot resize " << filename << " file to " << resize_bytes << " bytes\n"; return 1; }
	fclose(file);
	return 0;
}

int init_mmf(bi::file_mapping &fm, bi::mapped_region &mr, const wstring &wfilename, bi::mode_t mode, int64 resize_bytes) {
	if(resize_bytes > 0 && resize(wfilename, resize_bytes))
		return 1;
	// Current limitation of boost::interprocess, no wchar strings
	// std::codecvt
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
