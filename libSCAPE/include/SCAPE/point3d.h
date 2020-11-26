#pragma once

#include <iostream>
#include <cmath>

#include "point2d.h"


//template <bool, typename T>
//struct enable_if {
//};
//
//template <typename T>
//struct enable_if<true, T> {
//	typedef T result;
//};
//



//! 3D vector.
template <class T>
class point3d
{
public:

	static const unsigned int dim = 3;

	//! Value-Type of Templates components.
	typedef T ValueType; 
	typedef T value_type;
	typedef point2d<T> lowerDegreePoint;

// ****************************************************************************************

	point3d(ValueType v) {
		array[0] = array[1] = array[2] = v;
	}

// ****************************************************************************************

/*	point3d(int v) {
		array[0] = array[1] = array[2] = (ValueType)v;
	}
	*/

// ****************************************************************************************

	point3d() {
		array[0] = array[1] = array[2] = 0;
	}

// ****************************************************************************************

	point3d(T x, T y, T z) {
		array[0] = x;
		array[1] = y;
		array[2] = z;
	}

// ****************************************************************************************

	template <class U>
	point3d(const point3d<U>& other) {
		array[0] = (ValueType)other.array[0];
		array[1] = (ValueType)other.array[1];
		array[2] = (ValueType)other.array[2];
	}

// ****************************************************************************************

	point3d(const point3d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		array[2] = other.array[2];
	}

// ****************************************************************************************

	point3d(const T* other) {
		array[0] = other[0];
		array[1] = other[1];
		array[2] = other[2];
	}

// ****************************************************************************************

	inline bool operator<(const point3d& other) const {
		if ((x < other.x) && (y < other.y) && (z < other.z))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool lexicographicSmaller(const point3d& other) const {
		if (x < other.x) return true;
		if (x > other.x) return false;
		// x == other.x
		if (y < other.y) return true;
		if (y > other.y) return false;
		// also y == other.y
		if (z < other.z) return true;
		if (z > other.z) return false;

		// points are similar, thus NOT smaller
		return false;
	}

// ****************************************************************************************

	inline bool operator>(const point3d& other) const {
		if ((x > other.x) && (y > other.y) && (z > other.z))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator<=(const point3d& other) const {
		if ((x <= other.x) && (y <= other.y) && (z <= other.z))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator>=(const point3d& other) const {
		if ((x >= other.x) && (y >= other.y) && (z >= other.z))
			return true;

		return false;
	}

// ****************************************************************************************

	inline point3d<T> operator=(const point3d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		array[2] = other.array[2];
		return *this;
	}

// ****************************************************************************************

	//! Negieren
	inline point3d<T> operator-() const {
		return point3d<T>(-array[0], -array[1], -array[2]);
	}

// ****************************************************************************************

	inline point3d<T> operator+(const point3d& other) const {
		return point3d<T>(array[0]+other.array[0], array[1]+other.array[1], array[2]+other.array[2]);
	}

// ****************************************************************************************

	inline void operator+=(const point3d& other) {
		array[0] += other.array[0];
		array[1] += other.array[1];
		array[2] += other.array[2];
	}

// ****************************************************************************************

	inline void operator-=(const point3d& other) {
		array[0] -= other.array[0];
		array[1] -= other.array[1];
		array[2] -= other.array[2];
	}

// ****************************************************************************************

	inline void operator*=(T val) {
		array[0] *= val;
		array[1] *= val;
		array[2] *= val;
	}

// ****************************************************************************************

	inline void operator/=(T val) {

		T inv_val = ((T)1)/(val);

		array[0] *= inv_val;
		array[1] *= inv_val;
		array[2] *= inv_val;
	}

// ****************************************************************************************

	inline point3d<T> operator*(T val) const {
		return point3d<T>(array[0]*val, array[1]*val, array[2]*val);
	}

// ****************************************************************************************

	inline point3d<T> operator/(T val) const {
		return point3d<T>(array[0]/val, array[1]/val, array[2]/val);
	}

// ****************************************************************************************

	//! Vektor-/Kreuz-Produkt
	inline point3d<T> operator^(const point3d& other) const {
		return point3d<T>(array[1]*other.array[2] - array[2]*other.array[1], array[2]*other.array[0] - array[0]*other.array[2], array[0]*other.array[1] - array[1]*other.array[0]);
	}

// ****************************************************************************************

	//! Skalarmultiplikation
	inline T operator|(const point3d& other) const {
		return (array[0]*other.array[0] + array[1]*other.array[1] + array[2]*other.array[2]);
	}

// ****************************************************************************************

	//! Minus ;)
	inline point3d<T> operator-(const point3d& other) const {
		return point3d<T>(array[0]-other.array[0], array[1]-other.array[1], array[2]-other.array[2]);
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] == other[0]) && (this[1] == other[1])
	inline bool operator==(const point3d& other) const {
		if ((array[0] == other.array[0]) && (array[1] == other.array[1]) && (array[2] == other.array[2]))
			return true;

		return false;
	}

// ****************************************************************************************

	inline T squaredLength() const {
		return (array[0]*array[0] + array[1]*array[1] + array[2]*array[2]);
	}

// ****************************************************************************************

	inline T length() const {
		return sqrt(squaredLength());
	}

// ****************************************************************************************

	//! For marco compatibility
	inline T squaredDist(const point3d& other) const {
		return squareDist(other);
	}

// ****************************************************************************************

	inline T squareDist(const point3d& other) const {
		return ((array[0]-other.array[0])*(array[0]-other.array[0]) + (array[1]-other.array[1])*(array[1]-other.array[1]) + (array[2]-other.array[2])*(array[2]-other.array[2]));
	}

// ****************************************************************************************

	inline T dist(const point3d& other) const {
		return sqrt((((array[0]-other.array[0])*(array[0]-other.array[0]) + (array[1]-other.array[1])*(array[1]-other.array[1]) + (array[2]-other.array[2])*(array[2]-other.array[2]))));
	}

// ****************************************************************************************

	inline operator T*() {
		return array;
	}
	
// ****************************************************************************************

	inline operator const T*() const {
		return array;
	}
	
// ****************************************************************************************

	~point3d(void) {};

// ****************************************************************************************

	inline void print() const {
		std::cerr << "(" << array[0] << " " << array[1] << " " << array[2] << ")" << std::endl;
	}

//// ****************************************************************************************
//
//	const T& operator[](unsigned int i) const {
//		if (i < 3)
//			return array[i];
//		else {
//			std::cerr << "Out of bounds access to vec3" << std::endl;
//			exit(EXIT_FAILURE);
//		}
//	}
//
//// ****************************************************************************************
//
//	T& operator[](unsigned int i) {
//		if (i < 3)
//			return array[i];
//		else {
//			std::cerr << "Out of bounds access to vec3" << std::endl;
//			exit(EXIT_FAILURE);
//		}
//	}

// ****************************************************************************************

	inline point3d getNormalized() const {
		T val = length();
		return point3d<T>(array[0]/val, array[1]/val, array[2]/val);
	}

// ****************************************************************************************

	inline void normalize() {
		T val = length();
		array[0] /= val;
		array[1] /= val;
		array[2] /= val;
	}

// ****************************************************************************************

	union {
		struct {
			T x,y,z;          // standard names for components
		};
		struct {
			T r,g,b;          // standard names for components
		};
		T array[3];     // array access
	};
};


//! write a point3d to a stream
template <class T> inline std::ostream& operator<<(std::ostream& s, const point3d<T>& v)
{ return (s << v[0] << "/" << v[1] << "/" << v[2]);}

//! read a point3d from a stream
template <class T> inline std::istream& operator>>(std::istream& s, point3d<T>& v)
{ return (s >> v[0] >> v[1] >> v[2]); }


typedef point3d<double> vec3d;
typedef point3d<float> vec3f;
typedef point3d<int> vec3i;
typedef point3d<unsigned int> vec3ui;
typedef point3d<unsigned char> vec3uc;

//	
//template <> void point3d<double>::operator+=(const point3d<double>& other) {
//	x -= other.x;
//	y -= other.y;
//	z -= other.z;
//}



/*
Sicherheitskopie vor Umstellung von array[0] -> x, array[1] -> y, array[2] -> z.
*/