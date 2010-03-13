#pragma once
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <fstream>


#include <boost/cstdint.hpp>
typedef boost::int8_t int8;
typedef boost::int16_t int16;
typedef boost::int32_t int32;
typedef boost::uint8_t uint8;
typedef boost::uint16_t uint16;
typedef boost::uint32_t uint32;

typedef __int64 int64;
typedef unsigned __int64 uint64;

#include <boost/type_traits/make_signed.hpp>
//typedef size_t index;
//typedef boost::make_signed<index>::type signed_index;

#include <limits> // for maximum integer values
const double doubleNaN = std::numeric_limits<double>::quiet_NaN();
const float floatNaN = std::numeric_limits<float>::quiet_NaN();

#include <string>
#include <iostream>
#include <fstream>

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS // This is not enough sometimes, 
#pragma warning( disable : 4996)





enum { X, Y, Z };
//typedef float32 pt3float32[3];
//typedef float64 pt3float64[3];

// Types definitions
//template<int dim, typename t> struct Point;
template<typename t> struct Point3 {typedef t elem_type; t x,y,z;};
typedef Point3<float> float3;
typedef Point3<double> double3;


static inline istream &operator>>(istream &is, const char* str)
{
	string s(sizeof(str), ' ');
	is >> s;
	if(s != str)
		is.setstate(ios::badbit);
	return is;
}
static inline istream &operator>>(istream &is, const char ch)
{
	char ch_actual = is.get();
	return is;
}

template<typename t>
static inline ostream& operator <<(ostream &os, const Point3<t> &p) {
	os << '[' << p.x << ", " << p.y << ", " << p.z << ']';
	return os;
}
template<typename t>
static inline istream& operator>>(istream &is, Point3<t> &p) {
	is >> ws >> '[' >> p.x >> ws >> ',' >> p.y >> ws >>',' >> p.z >> ws >> ']';
	return is;
}


//Access function use: as_array(somexyzpoint)[2] is the same as somexyzpoint.z 
template<typename A>
inline A::elem_type * as_array(A &x) { return (A::elem_type*)&x; }

//Subscript index access
//template<int coord, typename point_type>
//inline point_type::elem_type & get(point_type &a) { return *((point_type::elem_type*)&a + coord); }


template<typename xyz>
inline xyz & operator+=(xyz &lhs, const xyz &rhs) {
	lhs.x += rhs.x;
	lhs.y += rhs.y;
	lhs.z += rhs.z;
	return lhs;
}
template<typename xyz>
inline xyz & operator-=(xyz &lhs, const xyz &rhs) {
	lhs.x -= rhs.x;
	lhs.y -= rhs.y;
	lhs.z -= rhs.z;
	return lhs;
}


template<typename xyz>
inline xyz operator-(const xyz &left, const xyz &right) {
	xyz res = {left.x - right.x, left.y - right.y, left.z - right.z};
	return res;
}

//template<typename float3, typename double3>
inline double3 operator+(const float3 &left, const double3 &right)
{
	double3 res = {left.x + right.x, left.y + right.y, left.z + right.z};
	return res;
}

template<typename xyz>
inline xyz & operator/=(xyz &lhs, const xyz::elem_type &rhs) {
	lhs.x /= rhs;
	lhs.y /= rhs;
	lhs.z /= rhs;
	return lhs;
}

inline float3 convert(const double3 &p) {
	float3 p1 = { (float)p.x, (float)p.y, (float)p.z};
	return p1;
};
inline double3 convert(const float3 &p) {
	double3 p1 = { p.x, p.y, p.z};
	return p1;
};

