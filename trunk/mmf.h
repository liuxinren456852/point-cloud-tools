#pragma once
namespace bi = boost::interprocess;
int init_mmf(bi::file_mapping &fm, bi::mapped_region &mr, const _TCHAR* file, bi::mode_t mode, int64 resize_bytes);

template<typename T>
class MMF : public boost::iterator_range<T*> {
	typedef boost::iterator_range<T*> base;
	bi::file_mapping fm;
	bi::mapped_region mr;
	wstring m_filename;
	//template<typename T> 
	T* pointer() const {return static_cast<T*>(mr.get_address());}
public:
 	int open(const _TCHAR* file, bi::mode_t mode, int64 resize_nelements) {
		m_filename = file;
		auto err = init_mmf(fm, mr, file, mode, resize_nelements*sizeof(value_type));
//		cerr << "File " << file << " will have " << mr.get_size()/sizeof(value_type) << "elements" << endl;
		if (!err)
			base::operator=(base(pointer(), pointer()+mr.get_size()/sizeof(value_type)));
		return err;
	}
	template<typename string_t> 
	int open(const string_t &f, bi::mode_t mode = bi::read_only, int64 resize = 0) { return open(f.c_str(), mode, resize); }
//	inline std::size_t get_size() const { return mr.get_size(); }
	void release() {  mr.swap(bi::mapped_region());  fm.swap(bi::file_mapping()); }
	int delete_file() {  release(); return _wremove(m_filename.c_str()); }
};
