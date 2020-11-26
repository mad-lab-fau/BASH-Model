#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "point1d.h"

//! 2D vector
template <class T>
class point2d
{
public:

	//! Value-Type of Templates components.
	typedef T ValueType;
	typedef T value_type;
	typedef point1d<T> lowerDegreePoint;

	static const unsigned int dim = 2;

// ****************************************************************************************

	point2d(T v) {
		array[0] = array[1] = v;
	}

// ****************************************************************************************

	point2d(const T* other) {
		array[0] = other[0];
		array[1] = other[1];
	}

// ****************************************************************************************

	point2d() {
		array[0] = array[1] = 0;
	}

// ****************************************************************************************

	point2d(T x, T y) {
		array[0] = x;
		array[1] = y;
	}

// ****************************************************************************************

	point2d(const point2d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
	}

// ****************************************************************************************

	point2d<T> operator=(const point2d& other) {
		array[0] = other.array[0];
		array[1] = other.array[1];
		return *this;
	}

// ****************************************************************************************

	point2d<T> operator-() {
		return point2d<T>(-array[0], -array[1]);
	}

// ****************************************************************************************

	point2d<T> operator+(const point2d& other) const {
		return point2d<T>(array[0]+other.array[0], array[1]+other.array[1]);
	}

// ****************************************************************************************

	void operator+=(const point2d& other) {
		array[0] += other.array[0];
		array[1] += other.array[1];
	}

// ****************************************************************************************

	void operator-=(const point2d& other) {
		array[0] -= other.array[0];
		array[1] -= other.array[1];
	}

// ****************************************************************************************

	void operator*=(T val) {
		array[0] *= val;
		array[1] *= val;
	}

// ****************************************************************************************

	void operator/=(T val) {
		array[0] /= val;
		array[1] /= val;
	}

// ****************************************************************************************

	point2d<T> operator*(T val) const {
		return point2d<T>(array[0]*val, array[1]*val);
	}

// ****************************************************************************************

	point2d<T> operator/(T val) const {
		return point2d<T>(array[0]/val, array[1]/val);
	}

// ****************************************************************************************

	inline point2d<T> operator-(const point2d& other) const {
		return point2d<T>(array[0]-other.array[0], array[1]-other.array[1]);
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] > other[0]) && (this[1] > other[1])
	bool operator>(const point2d& other)
	{
		if ((array[0] > other.array[0]) && (array[1] > other.array[1]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] >= other[0]) && (this[1] >= other[1])
	bool operator>=(const point2d& other)
	{
		if ((array[0] >= other.array[0]) && (array[1] >= other.array[1]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] <= other[0]) && (this[1] <= other[1])
	bool operator<=(const point2d& other) {
		if ((array[0] <= other.array[0]) && (array[1] <= other.array[1]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] == other[0]) && (this[1] == other[1])
	bool operator==(const point2d& other) {
		if ((array[0] == other.array[0]) && (array[1] == other.array[1]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] < other[0]) && (this[1] < other[1])
	bool operator<(const point2d& other) {
		if ((array[0] < other.array[0]) && (array[1] < other.array[1]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Skalarmultiplikation
	T operator|(const point2d& other) {
		return (array[0]*other.array[0] + array[1]*other.array[1]);
	}

// ****************************************************************************************

	T& operator[](unsigned int i) {
		if (i < 2)
			return array[i];
		else {
//			int j = i/0;
			std::cerr << "Out of bounds access to vec2" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

// ****************************************************************************************

	const T& operator[](unsigned int i) const {
		if (i < 2)
			return array[i];
		else {
			std::cerr << "Out of bounds access to vec2" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

// ****************************************************************************************

	~point2d(void) {};

// ****************************************************************************************

	inline bool operator<(const point2d& other) const {
		if ((x < other.x) && (y < other.y))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator>(const point2d& other) const {
		if ((x > other.x) && (y > other.y))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator<=(const point2d& other) const {
		if ((x <= other.x) && (y <= other.y))
			return true;

		return false;
	}

// ****************************************************************************************

	inline bool operator>=(const point2d& other) const {
		if ((x >= other.x) && (y >= other.y))
			return true;

		return false;
	}

// ****************************************************************************************

	T squaredLength() {
		return (array[0]*array[0] + array[1]*array[1]);
	}

// ****************************************************************************************

	//! Skalarmultiplikation
	inline T operator|(const point2d& other) const {
		return (array[0]*other.array[0] + array[1]*other.array[1]);
	}

// ****************************************************************************************

	T length() const {
		return sqrt(array[0]*array[0] + array[1]*array[1]);
	}
// ****************************************************************************************

	inline T squareDist(const point2d& other) const {
		return ((array[0]-other.array[0])*(array[0]-other.array[0]) + (array[1]-other.array[1])*(array[1]-other.array[1]));
	}

// ****************************************************************************************

	inline T dist(const point2d& other) const {
		return sqrt(((array[0]-other.array[0])*(array[0]-other.array[0]) + (array[1]-other.array[1])*(array[1]-other.array[1])));
	}

// ****************************************************************************************

	inline point2d getNormalized() const {
		T val = length();
		return point2d<T>(array[0]/val, array[1]/val);
	}

// ****************************************************************************************

	void normalize() {
		T val = length();
		array[0] /= val;
		array[1] /= val;
	}

// ****************************************************************************************

	void print() const {
		std::cerr << "(" << array[0] << " " << array[1] << ")" << std::endl;
	}

// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

	union {
		struct {
			T x,y;          // standard names for components
		};
		T array[2];     // array access
	};
};


//! write a point2d to a stream
template <class T> inline std::ostream& operator<<(std::ostream& s, const point2d<T>& v)
{ return (s << v[0] << " " << v[1]);}

//! read a point2d from a stream
template <class T> inline std::istream& operator>>(std::istream& s, point2d<T>& v)
{ return (s >> v[0] >> v[1]); }


typedef point2d<double> vec2d;
typedef point2d<float> vec2f;
typedef point2d<int> vec2i;
typedef point2d<unsigned int> vec2ui;
