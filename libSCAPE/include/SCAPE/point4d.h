#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "point3d.h"


//! 4D vector.
template <class T>
class point4d
{
public:

	//! Value-Type of Templates components.
	typedef T ValueType;
	typedef point3d<T> lowerDegreePoint;

	static const unsigned int dim = 4;

// ****************************************************************************************

	point4d(T v) {
		array[0] = array[1] = array[2] = array[3] = v;
	}

// ****************************************************************************************

	point4d() {
		array[0] = array[1] = array[2] = array[3] = 0;
	}

// ****************************************************************************************

	point4d(T x, T y, T z, T w) {
		array[0] = x;
		array[1] = y;
		array[2] = z;
		array[3] = w;
	}

// ****************************************************************************************

	point4d(const point3d<ValueType>& other, ValueType w) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		array[2] = other.array[2];
		array[3] = w;
	}

// ****************************************************************************************

	point4d(const point4d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		array[2] = other.array[2];
		array[3] = other.array[3];
	}

//// ****************************************************************************************
//
//	inline bool operator<(const point4d& other) const {
//		if ((x < other.x) && (y < other.y) && (z < other.z) && (w < other.w))
//			return true;
//
//		return false;
//	}

// ****************************************************************************************

	inline bool operator<(const point4d& other) const {
		if ((x < other.x) ||
			((x == other.x) && (y < other.y)) ||
			((x == other.x) && (y == other.y) && (z < other.z)) ||
			((x == other.x) && (y == other.y) && (z == other.z) && (w < other.w)))
			return true;

		return false;
	}
//


//// ****************************************************************************************
//
//	inline bool operator<(const point4d& other) const {
//		if ((x < other.x) && (y < other.y) && (z < other.z) && (w < other.w))
//			return true;
//
//		return false;
//	}

// ****************************************************************************************

	inline bool operator>(const point4d& other) const {
		if ((x > other.x) && (y > other.y) && (z > other.z) && (w > other.w))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator<=(const point4d& other) const {
		if ((x <= other.x) && (y <= other.y) && (z <= other.z) && (w <= other.w))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator>=(const point4d& other) const {
		if ((x >= other.x) && (y >= other.y) && (z >= other.z) && (w >= other.w))
			return true;

		return false;
	}

// ****************************************************************************************

	point4d<T> operator=(const point4d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		array[2] = other.array[2];
		array[3] = other.array[3];
		return *this;
	}

// ****************************************************************************************

	inline point4d<T> operator-() const {
		return point4d<T>(-array[0], -array[1], -array[2], -array[3]);
	}

// ****************************************************************************************

	inline point4d<T> operator+(const point4d& other) const {
		return point4d<T>(array[0]+other.array[0], array[1]+other.array[1], array[2]+other.array[2], array[3]+other.array[3]);
	}

// ****************************************************************************************

	inline void operator+=(const point4d& other) {
		array[0] += other.array[0];
		array[1] += other.array[1];
		array[2] += other.array[2];
		array[3] += other.array[3];
	}

// ****************************************************************************************

	inline void operator-=(const point4d& other) {
		array[0] -= other.array[0];
		array[1] -= other.array[1];
		array[2] -= other.array[2];
		array[3] -= other.array[3];
	}

// ****************************************************************************************

	inline void operator*=(T val) {
		array[0] *= val;
		array[1] *= val;
		array[2] *= val;
		array[3] *= val;
	}

// ****************************************************************************************

	inline void operator/=(T val) {
		array[0] /= val;
		array[1] /= val;
		array[2] /= val;
		array[3] /= val;
	}

// ****************************************************************************************

	inline point4d<T> operator*(T val) const {
		return point4d<T>(array[0]*val, array[1]*val, array[2]*val, array[3]*val);
	}

// ****************************************************************************************

	inline point4d<T> operator/(T val) const {
		return point4d<T>(array[0]/val, array[1]/val, array[2]/val, array[3]/val);
	}

// ****************************************************************************************

	//! Vektor-/Kreuz-Produkt
	inline point4d<T> operator^(const point4d& other) const {
		return point4d<T>(array[1]*other.array[2] - array[2]*other.array[1], array[2]*other.array[0] - array[0]*other.array[2], array[0]*other.array[1] - array[1]*other.array[0], T(1));
	}

// ****************************************************************************************

	//! Skalarmultiplikation
	inline T operator|(const point4d& other) const {
		return (array[0]*other.array[0] + array[1]*other.array[1] + array[2]*other.array[2] + array[3]*other.array[3]);
	}

// ****************************************************************************************

	inline point4d<T> operator-(const point4d& other) const {
		return point4d<T>(array[0]-other.array[0], array[1]-other.array[1], array[2]-other.array[2], array[3]-other.array[3]);
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] == other[0]) && (this[1] == other[1]) && ...
	inline bool operator==(const point4d& other) const {
		if ((array[0] == other.array[0]) && (array[1] == other.array[1]) && (array[2] == other.array[2]) && (array[3] == other.array[3]))
			return true;

		return false;
	}

// ****************************************************************************************

	inline T norm() const { return squaredLength(); }

// ****************************************************************************************

	inline T squaredLength() const {
		return (array[0]*array[0] + array[1]*array[1] + array[2]*array[2] + array[3]*array[3]);
	}

// ****************************************************************************************

	inline T length() const {
		return sqrt(squaredLength());
	}

// ****************************************************************************************

	inline T squareDist(const point4d& other) const {
		return (
			(array[0]-other.array[0])*(array[0]-other.array[0]) +
			(array[1]-other.array[1])*(array[1]-other.array[1]) +
			(array[2]-other.array[2])*(array[2]-other.array[2]) +
			(array[3]-other.array[3])*(array[3]-other.array[3])
			);
	}

// ****************************************************************************************

	inline T squaredDist(const point4d& other) const {
		return squareDist(other);
	}

// ****************************************************************************************

	inline T dist(const point4d& other) const {
		point4d<T> dist = *this - other;
		return dist.length();
	}

// ****************************************************************************************

	~point4d(void) {};

// ****************************************************************************************

	void print() const {
		std::cerr << "(" << array[0] << " " << array[1] << " " << array[2] << " " << array[3] << ")" << std::endl;
	}

// ****************************************************************************************

	inline const T& operator[](unsigned int i) const {

		return array[i];

		//if (i < 4)
		//	return array[i];
		//else {
		//	return array[i];
		//	std::cerr << "Out of bounds access to vec4" << std::endl;
		//	exit(EXIT_FAILURE);
		//}
	}

// ****************************************************************************************

	inline T& operator[](unsigned int i) {
		if (i < 4)
			return array[i];
		else {
			std::cerr << "Out of bounds access to vec4" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

// ****************************************************************************************

	inline void normalize() {
		T val = length();
		array[0] /= val;
		array[1] /= val;
		array[2] /= val;
		array[3] /= val;
	}

// ****************************************************************************************

	inline void dehomogenize() {
		array[0] /= array[3];
		array[1] /= array[3];
		array[2] /= array[3];
		array[3] /= array[3];
	}

// ****************************************************************************************

	inline bool isLinearDependent(const point4d& other) const {
		ValueType factor = x/other.x;

		if ((fabs(x/factor - other.x) + fabs(y/factor - other.y) + fabs(z/factor - other.z) + fabs(w/factor - other.w)) < 0.00001)
			return true;
		else
			return false;
	}

// ****************************************************************************************

	union {
		struct {
			T x,y,z,w;          // standard names for components
		};
		T array[4];     // array access
	};
};


//! write a point4d to a stream
template <class T> inline std::ostream& operator<<(std::ostream& s, const point4d<T>& v)
{ return (s << v[0] << "/" << v[1] << "/" << v[2] << "/" << v[3]);}

//! read a point4d from a stream
template <class T> inline std::istream& operator>>(std::istream& s, point4d<T>& v)
{ return (s >> v[0] >> v[1] >> v[2] >> v[3]); }


typedef point4d<double> vec4d;
typedef point4d<float> vec4f;
typedef point4d<int> vec4i;
