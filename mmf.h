#pragma once
namespace bi = boost::interprocess;
int init_mmf(bi::file_mapping &fm, bi::mapped_region &mr, const wstring &filename, bi::mode_t mode, int64 resize_bytes);
int resize(const wstring &filename, int64 resize_bytes);

template<typename T>
class MMF : public boost::iterator_range<T*> {
	typedef boost::iterator_range<T*> base;
	bi::file_mapping fm;
	bi::mapped_region mr;
	wstring m_filename;
	//template<typename T> 
	T* pointer() const {return static_cast<T*>(mr.get_address());}
public:
	MMF(const wstring filename) : m_filename(filename) {}
 	int open(bi::mode_t mode = bi::read_only, int64 resize_nelements = 0) {
		if(resize_nelements > 0)
			delete_file();
		auto err = init_mmf(fm, mr, m_filename, mode, resize_nelements*sizeof(value_type));
//		cerr << "File " << file << " will have " << mr.get_size()/sizeof(value_type) << "elements" << endl;
		if (!err)
			base::operator=(base(pointer(), pointer()+mr.get_size()/sizeof(value_type)));
		return err;
	}
	int rename_and_open() {
		const wstring newfn = m_filename + _T("_old");
		if(	_wrename(m_filename.c_str(), newfn.c_str())) {
			cerr << "Cannot rename the file\n";
			return 1;
		}
		m_filename = newfn;
		return open();
	}
//	template<typename string_t> 
//	int open(bi::mode_t mode = bi::read_only, int64 resize = 0) { return open(f.c_str(), mode, resize); }
//	inline std::size_t get_size() const { return mr.get_size(); }
	void release() {  mr.swap(bi::mapped_region());  fm.swap(bi::file_mapping()); }
	int delete_file() {  release(); return _wremove(m_filename.c_str()); }
	int resize_file(int64 newsize_nelements) {  release(); return resize(m_filename.c_str(), newsize_nelements*sizeof(value_type)); }
};
