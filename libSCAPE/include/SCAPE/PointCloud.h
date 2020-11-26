#pragma once

#include <vector>
#include <algorithm>
#include <limits>

template<class T>
struct sort_pred {
	bool operator () (const T& left, const T& right) {
		//return left.z < right.z;
		return left[T::dim-1] < right[T::dim-1];
	}
};



//! This class represents a pointcloud or unstructured pointset, i.e. it holds nothing but a vector of 2D/3D points
template<class T>
class PointCloud {

public:

	//! Type of a Point, for example a 3D Vector of floats.
	typedef T Point;

	//! Value-Type of Templates components (i.e. float or double).
	typedef typename T::ValueType ValueType;

	//! Dimension of the pointcloud.
	static const UInt dim = T::dim;

	PointCloud() {
		maxPoint = Point(-std::numeric_limits<ValueType>::max());
		minPoint = Point(std::numeric_limits<ValueType>::max());
	}
	
	virtual ~PointCloud() {}

	//! Recomputes all vertices as pt_new = (pt+move)*scale;
	inline void resize(ValueType scale, Point move) {
		for (int i1 = 0; i1 < getNumPts(); i1++) {
			m_pts[i1] = (m_pts[i1]+move)*scale;
		}
	}

	void recomputeMinMax() {
		maxPoint = Point(-std::numeric_limits<ValueType>::max());
		minPoint = Point(std::numeric_limits<ValueType>::max());

		for (int i1 = 0; i1 < getNumPts(); i1++) {
			for (UInt i = 0; i < dim; ++i) {
				minPoint[i] = std::min(minPoint[i], m_pts[i1][i]);
				maxPoint[i] = std::max(maxPoint[i], m_pts[i1][i]);
			}
		}
	}
	

	//! Add a point to the pointcloud
	virtual void insertPoint(const T& pt) {
		for (UInt i = 0; i < dim; ++i) {
			minPoint[i] = std::min(minPoint[i], pt[i]);
			maxPoint[i] = std::max(maxPoint[i], pt[i]);
		}
		m_pts.push_back(pt);
	};

	//! Gets the number of points in the pointcloud
	int getNumPts() const {
		return (int)m_pts.size();
	};

	//! Returns an array of the points
	const Point* getPoints() const {
		return &m_pts[0];
	};

	//! Returns an array of the points
	Point* getPoints() {
		return &m_pts[0];
	};

	//! Returns the point with index 'i'
	inline Point getPoint(int i) {
		return m_pts[i];
	};

	//! Lower left front point of the extend.
	Point getMin() {
		return minPoint;
	};

	//! Upper right back point of the extend.
	Point getMax() {
		return maxPoint;
	};

	//! Prints all points.
	void print() {
		for (UInt i1 = 0; i1 < m_pts.size(); i1++) {
			m_pts[i1].print();
		}
	}

	void normalizeConservative(ValueType& scale_factor, Point& move) {
		scale_factor = 0;
		for (UInt i1 = 0; i1 < dim; i1++) {
			scale_factor = std::max((maxPoint[i1]- minPoint[i1]), scale_factor);
		}

		scale_factor *= (ValueType)1.05;

		Point min_tmp = minPoint/scale_factor;
		Point max_tmp = maxPoint/scale_factor;
		move = ((max_tmp+min_tmp)/(ValueType)2)-Point((ValueType)0.5);

		for (UInt i1 = 0; i1 < m_pts.size(); i1++) {
			m_pts[i1] = (m_pts[i1]/scale_factor) - move;
			//m_pts[i1].print();
		}

		recomputeMinMax();
		//std::cerr << minPoint << "  " << maxPoint << std::endl;
	}

	void normalize(ValueType& scale_factor, Point& move) {
		scale_factor = 0;
		for (UInt i1 = 0; i1 < dim; i1++) {
			scale_factor = std::max((maxPoint[i1] - minPoint[i1]), scale_factor);
		}

		Point min_tmp = minPoint/scale_factor;
		Point max_tmp = maxPoint/scale_factor;
		move = ((max_tmp+min_tmp)/(ValueType)2)-Point((ValueType)0.5);

		for (UInt i1 = 0; i1 < m_pts.size(); i1++) {
			m_pts[i1] = (m_pts[i1]/scale_factor) - move;
			//m_pts[i1].print();
		}

		recomputeMinMax();
		//std::cerr << minPoint << "  " << maxPoint << std::endl;
	}

	//! Sorts the point in z-dir.
	void sort_Z() {
		std::sort(m_pts.begin(), m_pts.end(), sort_pred<T>());
	}

	virtual void null() {};

	//! Storage for the points of the pointset.
	std::vector<T> m_pts;

	//! Extend of the pointcloud, min[0] denotes the smallest x value, max[1] the largest y value...
	T minPoint, maxPoint;

protected:


};



//****************************************************************************************************


//! This class represents a pointcloud with normals
template<class T>
class PointCloudNormals : public PointCloud<T> {

public:
	// FIX: https://stackoverflow.com/questions/1567730/inheritance-and-templates-in-c-why-are-inherited-members-invisible
	using RigidTransform<T>::Point;
	using RigidTransform<T>::ValueType;
	using RigidTransform<T>::minPoint;
	using RigidTransform<T>::maxPoint;
	using RigidTransform<T>::m_pts;

	PointCloudNormals(): PointCloud() {
	}

	//! Destruktor, frees some Memory
	PointCloudNormals(const PointCloudNormals& other) {
		minPoint = other.minPoint;
		maxPoint = other.maxPoint;
		m_pts = other.m_pts;
		m_normals = other.m_normals;
	}

	//! Destruktor, frees some Memory
	~PointCloudNormals() {
		m_pts.clear();
		m_normals.clear();
	}

	//! Add a point to the pointcloud
	void insertPoint(const T& pt, const T& normal) {
		PointCloud<T>::insertPoint(pt);
		m_normals.push_back(normal);
	}

	//! Returns an array of the points
	const T* getNormals() const{
		return &m_normals[0];
	};
	
	//! Returns an array of the points
	T* getNormals() {
		return &m_normals[0];
	};

	//! Storage for the normals of the pointset.
	std::vector<T> m_normals;

private:

};
