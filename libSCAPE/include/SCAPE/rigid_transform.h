#pragma once

#include "MatrixnD.h"
//#include "SimpleMesh.h"

template <class T>
class RigidTransform {
public:

	//! Type of a Vector.
	typedef T Point;

	//! floating point type.
	typedef typename Point::ValueType ValueType;

	//! Dimension nxn.
	static const int dim = Point::dim;

	//! Typedef for a matrix of the current type.
	typedef SquareMatrixND<T> Matrix;

	//! Constructor, sets everything to identity.
	RigidTransform() {
		rotate.setToIdentity();
		translate = Point(ValueType(0));
	}

	//! Constructs a rigid transformation from a rotation matrix and a translation.
	RigidTransform(const Point& t, const Matrix& r) {
		translate = t;
		rotate = r;
	}

	//! Returns the inverse rigid transformation.
	RigidTransform getInverseTransform() const {
		RigidTransform<T> inv;

		inv.rotate = rotate;
		inv.rotate.transpose();

		inv.translate = inv.rotate.vecTrans(-translate);

		return inv;
	}

	//! Transforms a normal (i.e. skips the translation part).
	inline Point normTrans(const Point& t) const {
		return rotate.vecTrans(t);
	}

	//! Transforms a vector.
	inline Point vecTrans(const Point& t) const {
		Point t2 = rotate.vecTrans(t);
		t2 += translate;
		return t2;
	}

	//! Transforms all elements in the vector
	void vecArrayTrans(std::vector<Point>& vecs) const {
		for (unsigned int i1 = 0; i1 < vecs.size(); i1++) vecs[i1] = vecTrans(vecs[i1]);
	}


	//! Appends a rigid transformation to the current rigid transformation.
	void appendRigidTransformation(const RigidTransform<T>& other) {
		rotate = other.rotate*rotate;
		translate = other.rotate.vecTrans(translate)+other.translate;
	}

	//! Prints out rotation-matrix and translation-vector.
	void print() const {
		std::cerr << "Rotation: " << std::endl;
		rotate.print();
		std::cerr << "Translation: " << std::endl;
		translate.print();
	}

	Matrix getRot() const {
		return rotate;
	}

	Point getTrans() const {
		return translate;
	}

protected:

public:
	//! Translation.
	Point translate;
	//! Rotation.
	Matrix rotate;

};

template <class T>
class RigidTransformWithScale : public RigidTransform<T> {

public:
	// FIX: https://stackoverflow.com/questions/1567730/inheritance-and-templates-in-c-why-are-inherited-members-invisible
	using RigidTransform<T>::Point;
	using RigidTransform<T>::ValueType;
	using RigidTransform<T>::Matrix;
	using RigidTransform<T>::translate;
	using RigidTransform<T>::rotate;

	//! Constructor, sets everything to identity.
	RigidTransformWithScale() : RigidTransform()  {
		scale = 0;
	}

	//! Constructs a rigid transformation from a rotation matrix and a translation.
	RigidTransformWithScale(const Point& t, const Matrix& r, const ValueType& s) : RigidTransform<T>(t, r) {
		scale = s;
		//std::cerr << "Scale: " << scale << std::endl;
	}

	//! Prints out rotation-matrix and translation-vector.
	void print() {
		std::cerr << "Rotation: " << std::endl;
		rotate.print();
		std::cerr << "Translation: " << std::endl;
		translate.print();
		std::cerr << "Scale: " << scale << std::endl;
	}

	//! Returns the inverse rigid transformation.
	RigidTransformWithScale getInverseTransform() const {
		RigidTransformWithScale<T> inv;

		inv.rotate = rotate;
		inv.rotate.transpose();

		inv.translate = inv.rotate.vecTrans(-translate)/scale;

		inv.scale = (ValueType)1.0 / scale;

		return inv;
	}

	inline void setScale(ValueType scale_) {
		scale = scale_;
	}

	//! Transforms a normal (i.e. skips the translation part).
	inline Point normTrans(const Point& t) const {
		return rotate.vecTrans(t);
	}

	//! Transforms a vector.
	inline Point vecTrans(const Point& t) {
		Point t2 = rotate.vecTrans(t*scale);
		return t2+translate;
	}

	//! Appends a rigid transformation to the current rigid transformation.
	void appendRigidTransformation(const RigidTransformWithScale<T>& other) {
		scale *= other.scale;
		rotate = other.rotate*rotate;
		translate = other.rotate.vecTrans(translate)+other.translate;
	}

	//! Scaling
	ValueType scale;

private:

};