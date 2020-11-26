#pragma once

// Disables the 'nameless struct/union' warning in VC.
#pragma warning( disable : 4201 )


#include <iostream>
#include <cmath>
#include <cstdlib>

//! 1D vector (I know it's a joke, but we need it for compatibility reasons)
template <class T>
class point1d
{
public:

	//! Value-Type of Templates components.
	typedef T ValueType;

	static const unsigned int dim = 1;

// ****************************************************************************************

	point1d(T v) {
		array[0] = v;
	}

// ****************************************************************************************

	point1d() {
		array[0] = 0;
	}

// ****************************************************************************************

	point1d(const point1d& other) {
		array[0] = other.array[0];
	}

// ****************************************************************************************

	point1d<T> operator=(const point1d& other) {
		array[0] = other.array[0];
		return *this;
	}

// ****************************************************************************************

	point1d<T> operator-() {
		return point1d<T>(-array[0]);
	}

// ****************************************************************************************

	point1d<T> operator+(const point1d& other) {
		return point1d<T>(array[0]+other.array[0]);
	}

// ****************************************************************************************

	void operator+=(const point1d& other) {
		array[0] += other.array[0];
	}

// ****************************************************************************************

	void operator-=(const point1d& other) {
		array[0] -= other.array[0];
	}

// ****************************************************************************************

	void operator*=(T val) {
		array[0] *= val;
	}

// ****************************************************************************************

	void operator/=(T val) {
		array[0] /= val;
	}

// ****************************************************************************************

	point1d<T> operator*(T val) {
		return point1d<T>(array[0]*val);
	}

// ****************************************************************************************

	point1d<T> operator/(T val) {
		return point1d<T>(array[0]/val);
	}

// ****************************************************************************************

	point1d<T> operator-(const point1d& other) {
		return point1d<T>(array[0]-other.array[0]);
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] > other[0]).
	bool operator>(const point1d& other)
	{
		if ((array[0] > other.array[0]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] >= other[0]).
	bool operator>=(const point1d& other)
	{
		if ((array[0] >= other.array[0]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] <= other[0]).
	bool operator<=(const point1d& other) {
		if ((array[0] <= other.array[0]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] == other[0]).
	bool operator==(const point1d& other) {
		if ((array[0] == other.array[0]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Liefert nur dann true, wenn (this[0] < other[0]) && (this[1] < other[1])
	bool operator<(const point1d& other) {
		if ((array[0] < other.array[0]))
			return true;

		return false;
	}

// ****************************************************************************************

	//! Skalarmultiplikation
	T operator|(const point1d& other) {
		return (array[0]*other.array[0]);
	}

// ****************************************************************************************

	T& operator[](unsigned int i) {
		if (i < 1)
			return array[i];
		else {
			std::cerr << "Out of bounds access to vec1" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

// ****************************************************************************************

	const T& operator[](unsigned int i) const {
		if (i < 1)
			return array[i];
		else {
			std::cerr << "Out of bounds access to vec1" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

// ****************************************************************************************

	~point1d(void) {};

// ****************************************************************************************

	T squaredLength() {
		return (array[0]*array[0]);
	}

// ****************************************************************************************

	T length() {
		return array[0];
	}

// ****************************************************************************************

	T dist(const point1d& other) {
		point1d<T> dist = *this - other;
		return dist.length();
	}

// ****************************************************************************************

	void normalize() {
		array[0] = 1;
	}

// ****************************************************************************************

	void print() const {
		std::cerr << "(" << array[0] << ")" << std::endl;
	}

// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

	union {
		struct {
			T x;          // standard names for components
		};
		T array[1];     // array access
	};
};


//! write a point1d to a stream
template <class T> inline std::ostream& operator<<(std::ostream& s, const point1d<T>& v)
{ return (s << v[0]);}

//! read a point1d from a stream
template <class T> inline std::istream& operator>>(std::istream& s, point1d<T>& v)
{ return (s >> v[0]); }


typedef point1d<double> vec1d;
typedef point1d<float> vec1f;
typedef point1d<int> vec1i;
