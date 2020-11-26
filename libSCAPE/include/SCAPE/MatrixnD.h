#pragma once


#include "point2d.h"
#include "point3d.h"
#include "nr_templates.h" // for eigenvectors/eigenvalues (jacobi)
#include <algorithm>
#include <iomanip>
#include <vector>


template<class T>
class FloatingPointIntPair {
public:
	typedef T ValueType;

	FloatingPointIntPair(int i_, ValueType l_) : index(i_), lambda(l_) {}

	inline bool operator<(const FloatingPointIntPair& other) const {
		if (lambda > other.lambda)
			return true;

		return false;
	}

	int index;
	ValueType lambda;
};

template<class T> class SquareMatrixND;

template <class T>
inline T determinant3x3(const SquareMatrixND< point3d<T> >& mat) {
	T res = 0;

	res += mat.m[0]*mat.m[4]*mat.m[8];
	res += mat.m[1]*mat.m[5]*mat.m[6];
	res += mat.m[2]*mat.m[3]*mat.m[7];
	res -= mat.m[2]*mat.m[4]*mat.m[6];
	res -= mat.m[1]*mat.m[3]*mat.m[8];
	res -= mat.m[0]*mat.m[5]*mat.m[7];
	return res;
}

// ****************************************************************************************


//! A n-dimensional Square-Matrix
/*!
	Ordering (for operator()) and set/get:<br>
	0/0 0/1 0/2 ... 0/n<br>
	1/0 1/1 1/2 ... 1/n<br>
	...<br>
	n/0 n/1 n/2 ... n/n<br>
*/
template<class T>
class SquareMatrixND {

public:

	typedef unsigned int UInt;
	typedef T Point;
	typedef typename Point::ValueType ValueType;
	static const int dim = Point::dim;
	static const int numFields = dim*dim;

// ****************************************************************************************

	//! Constructs an empty SquareMatrix.
	SquareMatrixND() {
		for (UInt i1 = 0; i1 < numFields; i1++) m[i1] = 0;
	}

// ****************************************************************************************

	//! Constructs matrix from the given array (assumed to be of size dim*dim!).
	SquareMatrixND(const ValueType* array) {
		for (UInt i1 = 0; i1 < numFields; i1++) m[i1] = array[i1];
	}

// ****************************************************************************************

	//! Constructs an const  SquareMatrix.
	SquareMatrixND(ValueType val) {
		for (UInt i1 = 0; i1 < numFields; i1++) m[i1] = val;
	}

// ****************************************************************************************

	//! Copy Constructor.
	SquareMatrixND(const SquareMatrixND<Point>& other) {
		for (UInt i1 = 0; i1 < numFields; i1++) {
			m[i1] = other.m[i1];
		}
	}

// ****************************************************************************************

	//! Constructs a SquareMatrix from the Tensor-Product of 'vec'
	SquareMatrixND(const Point& vec) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = i1; i2 < dim; i2++) {
				m[i1*dim+i2] = vec[i1]*vec[i2];
			}
		}
		symmetrize();
	}

// ****************************************************************************************

	//! Returns true if all entries are zero.
	inline bool isZero() const {
		for (UInt i1 = 0; i1 < numFields; i1++) if (m[i1] != 0) return false;
		return true;
	}

// ****************************************************************************************

	//! Sets all elements to zero.
	inline void zeroOut() {
		for (UInt i1 = 0; i1 < numFields; i1++) m[i1] = 0;
	}

// ****************************************************************************************

	//! Overwrites current matrix with the identity matrix.
	inline void setToIdentity() {
		for (UInt i1 = 0; i1 < numFields; i1++) m[i1] = 0;
		for (UInt i1 = 0; i1 < dim; i1++) (*this)(i1,i1) = 1;
	}


// ****************************************************************************************

	//! Returns the rotation vector of the current matrix. Inverse function to 'setToRotationMatrixNew'.
	inline Point getRotVec() const {

		Point r((*this)(1,2)-(*this)(2,1), (*this)(2,0)-(*this)(0,2), (*this)(0,1)-(*this)(1,0));
		ValueType len = r.length();
		if (len == 0) return r;
		r *= (atan2(len, (*this)(0,0)+(*this)(1,1)+(*this)(2,2)-1)/len);
		return r;

	}

// ****************************************************************************************

	//! Creates a rotation-matrix from an axis and an angle (Works for 3x3 matrices only!). Inverse function to 'getRotVec'.
	inline void setToRotationMatrixNew(Point c) {

		if (c.squaredLength() == 0) {
			setToIdentity();
			return;
		}

		ValueType theta = c.length();
		c /= theta;

		SquareMatrixND<Point> R;
		R(0,1) = c.z;
		R(0,2) = -c.y;
		R(1,0) = -c.z;
		R(1,2) = c.x;
		R(2,0) = c.y;
		R(2,1) = -c.x;

		(*this) = (R * sin(theta) + (R*R)*(1-cos(theta)));
		(*this)(0,0) += 1;
		(*this)(1,1) += 1;
		(*this)(2,2) += 1;

	}

// ****************************************************************************************

	//! Creates (caution, gives the transpose=inverted result) a rotation-matrix from an axis and an angle (Works for 3x3 matrices only!).
	inline void setToRotationMatrix_deprecated(Point c, ValueType a) {
		c.normalize();
		ValueType sa = sin(a);
		ValueType ca = cos(a);
		ValueType ca2 = 1-ca;
		(*this)(0,0) = ca			+	c.x*c.x*ca2;
		(*this)(0,1) = c.x*c.y*ca2	-	c.z*sa;
		(*this)(0,2) = c.x*c.z*ca2	+	c.y*sa;
		(*this)(1,0) = c.y*c.x*ca2	+	c.z*sa;
		(*this)(1,1) = ca			+	c.y*c.y*ca2;
		(*this)(1,2) = c.y*c.z*ca2	-	c.x*sa;
		(*this)(2,0) = c.z*c.x*ca2	-	c.y*sa;
		(*this)(2,1) = c.z*c.y*ca2	+	c.x*sa;
		(*this)(2,2) = ca			+	c.z*c.z*ca2;
	}

// ****************************************************************************************

	//! Copy the values from the upper left matrix half to the lower right
	inline void symmetrize() {
		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = i1+1; i2 < dim; i2++) {
				m[i2*dim+i1] = m[i1*dim+i2];
			}
		}
	}

// ****************************************************************************************

	//! Adds the Tensor-Product of 'vec1' x 'vec2' to the matrix.
	inline void addFromTensorProduct(const Point& vec1, const Point& vec2) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				m[i1*dim+i2] += vec1[i1]*vec2[i2];
			}
		}
	}

// ****************************************************************************************

	//! Adds the Tensor-Product of 'vec' to the UPPER-LEFT half of the matrix, does not affect the lower half!!! <b>call symmetrize when done!!</b>.
	inline void addFromTensorProduct_NO_SYMMETRIZE(const Point& vec) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = i1; i2 < dim; i2++) {
				m[i1*dim+i2] += vec[i1]*vec[i2];
			}
		}
	}

// ****************************************************************************************

	//! Sets the 'i'th column of the matrix to vec.
	inline void setColI(const Point& vec, UInt i) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			m[i1*dim+i] = vec[i1];
		}
	}

// ****************************************************************************************

	//! sets the values in the 'i'th row of the matrix to 'vec'
	inline void setRowI(const Point& vec, UInt i) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			m[i*dim+i1] = vec[i1];
		}
	}

// ****************************************************************************************

	//! Returns the diagonal of the matrix
	inline Point diag() const {
		Point ret;
		for (UInt i1 = 0; i1 < dim; i1++) ret[i1] = m[i1*dim+i1];
		return ret;
	}

// ****************************************************************************************

	//! Access element 'i','j'.
	inline const ValueType& operator() (const int i, const int j) const {
		return m[i*dim+j];	
	};

// ****************************************************************************************

	//! Access element 'i','j'.
	inline ValueType& operator() (const int i, const int j) {
		return m[i*dim+j];	
	};

// ****************************************************************************************

	//! Multiplies all entries by the factor 'fac'.
	void operator*=(const ValueType fac) {
		for (UInt i1 = 0; i1 < numFields; i1++) {
			m[i1]*=fac;
		}
	}

// ****************************************************************************************

	//! Divides all entries by the factor 'fac'.
	void operator/=(const ValueType fac) {
		for (UInt i1 = 0; i1 < numFields; i1++) {
			m[i1]/=fac;
		}
	}

// ****************************************************************************************

	//! returns (*this)-'other'.
	SquareMatrixND<Point> operator-(const SquareMatrixND<Point>& other) const {
		SquareMatrixND<Point> ret;
		for (UInt i1 = 0; i1 < numFields; i1++) {
			ret.m[i1] = m[i1]-other.m[i1];
		}
		return ret;
	}

// ****************************************************************************************

	//! returns (*this)+'other'.
	SquareMatrixND<Point> operator+(const SquareMatrixND<Point>& other) const {
		SquareMatrixND<Point> ret;
		for (UInt i1 = 0; i1 < numFields; i1++) {
			ret.m[i1] = m[i1]+other.m[i1];
		}
		return ret;
	}

// ****************************************************************************************

	//! returns (*this)*'fac'.
	SquareMatrixND<Point> operator*(const ValueType fac) const {
		SquareMatrixND<Point> ret;
		for (UInt i1 = 0; i1 < numFields; i1++) {
			ret.m[i1] = m[i1]*fac;
		}
		return ret;
	}

// ****************************************************************************************

	//! Adds 'other' to the current matrix.
	void operator+=(const SquareMatrixND<Point>& other) {
		for (UInt i1 = 0; i1 < numFields; i1++) {
			m[i1] += other.m[i1];
		}
	}

// ****************************************************************************************

	//! Matrixmultiplication: returns (*this)*'other'.
	SquareMatrixND<Point> operator*(const SquareMatrixND<Point>& other) const {

		SquareMatrixND<Point> ret;

		for(UInt i1 = 0; i1 < dim; i1++) {
			for(UInt i2 = 0; i2 < dim; i2++)	{
				ret(i1,i2) = 0.0;
				for(UInt i3 = 0; i3 < dim; i3++)
					ret(i1,i2) += (*this)(i1,i3) * other(i3,i2);
				
			}
		}
		return ret;
	}
	
// ****************************************************************************************

	//! writes the entries of row 'i' to row 'j' and vice versa.
	inline void changeRows(UInt i, UInt j) {
		Point row_i = getRowI(i);
		Point row_j = getRowI(j);
		setRowI(row_j, i);
		setRowI(row_i, j);
	}

// ****************************************************************************************

	//! returns the 'i'th row of the matrix
	inline void flipSignsInRowI(UInt i) {
		for (UInt i1 = 0; i1 < dim; i1++) {
			m[i*dim+i1] = -m[i*dim+i1];
		}
	}

// ****************************************************************************************

	//! returns the 'i'th row of the matrix
	inline Point getRowI(UInt i) const {
		Point row;
		for (UInt i1 = 0; i1 < dim; i1++) {
			row[i1] = m[i*dim+i1];
		}
		return row;
	}

// ****************************************************************************************

	//! returns the 'i'th column of the matrix
	inline Point getColI(UInt i) const {
		Point row;
		for (UInt i1 = 0; i1 < dim; i1++) {
			row[i1] = m[i1*dim+i];
		}
		return row;
	}

// ****************************************************************************************

	//! Prints matrix to stderr.
	void print(int precision = 5) const {
		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				std::cerr << std::fixed << std::setprecision(precision) << m[i1*dim+i2] << " ";
			}
			std::cerr << std::endl;
		}
	}

// ****************************************************************************************

	//! Calculates the eigenvalues and eigenvectors of the matrix. Note that the eigensystem is not necessarily right handed!
	Point calcEValuesAndVectorsCORRECT(SquareMatrixND<Point>& ESystem) {

		// Use jacobi's method:
		// Build dim x dim matrix NR-style:
		float** CV = new float*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++)
			CV[i1] = new float[dim+1];
		float lambda[dim+1];
		float** v = new float*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++)
			v[i1] = new float[dim+1];


		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				CV[i+1][j+1] = (float)(*this)(i,j);
			}
		}

		int num_of_required_jabobi_rotations;

		if (!jacobi(CV, dim, lambda, v, &num_of_required_jabobi_rotations)) {
			std::cerr << "ERROR: could not compute Eigenvectors in SquareMatrixND::calcEValuesAndVectors(...)" << std::endl;
			print();
			std::cerr << "EXIT" << std::endl;
			exit(EXIT_FAILURE);
		}

		Point eigenvalues;

		// sort EV's correctly:
		std::vector< FloatingPointIntPair<ValueType> > llps;
		for (UInt i1 = 1; i1 <= dim; i1++) {
			llps.push_back( FloatingPointIntPair<ValueType>(i1, lambda[i1]) );
		}
		std::sort(llps.begin(), llps.end());

		//for (UInt i1 = 0; i1 < dim; i1++) {
		//	std::cerr << llps[i1].lambda << std::endl;
		//}

		for (UInt i1 = 0; i1 < dim; i1++) {
			// Build vector:
			Point p;
			for (UInt i2 = 0; i2 < dim; i2++) {
				p[i2] = v[i2+1][llps[i1].index];
			}
			eigenvalues[i1] = llps[i1].lambda;
			ESystem.setColI(p, i1);
			//checkIsEigenvector(p);
		}

		// clean up
		for(int i1 = 0; i1 < dim+1; i1++) delete[] CV[i1];
		delete[] CV;

		for(int i1 = 0; i1 < dim+1; i1++) delete[] v[i1];
		delete[] v;

		return eigenvalues;
	}

// ****************************************************************************************

	//! Calculates the eigenvalues and eigenvectors of the matrix. Note that the eigensystem is not necessarily right handed!
	Point calcEValuesAndVectors(SquareMatrixND<Point>& ESystem) {

		// Use jacobi's method:
		// Build dim x dim matrix NR-style:
		float** CV = new float*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++)
			CV[i1] = new float[dim+1];
		float lambda[dim+1];
		float** v = new float*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++)
			v[i1] = new float[dim+1];


		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				CV[i+1][j+1] = (float)(*this)(i,j);
			}
		}

		int num_of_required_jabobi_rotations;

		if (!jacobi(CV, dim, lambda, v, &num_of_required_jabobi_rotations)) {
			std::cerr << "ERROR: could not compute Eigenvectors in SquareMatrixND::calcEValuesAndVectors(...)" << std::endl;
			print();
			std::cerr << "EXIT" << std::endl;
			exit(EXIT_FAILURE);
		}

		Point eigenvalues;

		//std::cerr << std::endl << std::endl;
		//for (int i = 1; i <= dim; i++) {
		//	for (int j = 1; j <= dim; j++) {
		//		std::cerr << v[i][j] << "   ";
		//	}
		//	std::cerr << "  + " << lambda[i] << std::endl;
		//}
		//std::cerr << std::endl << std::endl;

		// sort EV's correctly:
		std::vector< FloatingPointIntPair<ValueType> > llps;
		for (UInt i1 = 1; i1 <= dim; i1++) {
			llps.push_back( FloatingPointIntPair<ValueType>(i1, lambda[i1]) );
		}
		std::sort(llps.begin(), llps.end());

		//for (UInt i1 = 0; i1 < dim; i1++) {
		//	std::cerr << llps[i1].lambda << std::endl;
		//}

		for (UInt i1 = 0; i1 < dim; i1++) {
			// Build vector:
			Point p;
			for (UInt i2 = 0; i2 < dim; i2++) {
				p[i2] = v[i2+1][llps[i1].index];
			}
			eigenvalues[i1] = (llps[i1].lambda != 0) ? (ValueType)1.0/llps[i1].lambda : 1000000;
			//eigenvalues[i1] = llps[i1].lambda;
			ESystem.setColI(p, i1);
			//checkIsEigenvector(p);
		}

		// clean up
		for(int i1 = 0; i1 < dim+1; i1++) delete[] CV[i1];
		delete[] CV;

		for(int i1 = 0; i1 < dim+1; i1++) delete[] v[i1];
		delete[] v;

		return eigenvalues;
	}

// ****************************************************************************************

	void checkIsEigenvector(const Point& vec) const {
		Point vec_trans = vecTrans(vec);
		vec_trans.normalize();
		std::cerr << "cIE: " << (vec|vec_trans) << std::endl;
	}

// ****************************************************************************************

	//! Returns the inverted matrix.
	inline SquareMatrixND<Point> getInverted(ValueType eps = (ValueType)0.00000001) const {

		//std::cerr << "Matrix inversion is evil" << std::endl;

		SquareMatrixND<Point> U, V;
		Point sigma;
		SVD_decomp(U, sigma, V);

		// build pseudo-inverse:
		SquareMatrixND<Point> S;
		for (UInt i1 = 0; i1 < dim; i1++) {
			if (sigma[i1] < eps) S(i1,i1) = 0;
			else S(i1,i1) = ((ValueType)1.0)/sigma[i1];
		}

		U.transpose();

		return V*(S*U);

	}

// ****************************************************************************************

	//! Returns the best fitting orthonormal matrix 'ortho' such that ||'this'-'ortho'||fröbenius == min
	SquareMatrixND<Point> getBestOrtho() const {

		SquareMatrixND<Point> U, V;
		Point sigma;

		SVD_decomp(U, sigma, V);
		V.transpose();
		SquareMatrixND<Point> ortho = U*V;

		double det = ortho.getDeterminant();
		if (det < 0) { // SEEMS NOT NECESSARY
			int sm_id = 0;
			double smallest = sigma[sm_id];
			for (UInt dd = 1; dd < 3; dd++) {
				if (sigma[dd] < smallest) {
					smallest = sigma[dd];
					sm_id = dd;
				}
			}
			// flip sign of entries in colums 'sm_id' in 'U'
			U.m[sm_id + 0] *= -1;
			U.m[sm_id + 3] *= -1;
			U.m[sm_id + 6] *= -1;
			ortho = U*V;
		}

		
		//// compute / check fröbenius norm
		//SquareMatrixND<Point> diff = (*this)-ortho;
		//diff.SVD_decomp(U, sigma, V);
		//ValueType s_l = sigma[0];
		//if (sigma[1] > s_l) s_l = sigma[1];
		//if (sigma[2] > s_l) s_l = sigma[2];
		//std::cerr << s_l << std::endl;

		return ortho;
		
	}

// ****************************************************************************************

	//! Returns the transposed matrix.
	inline SquareMatrixND<Point> getTransposed() const {

		SquareMatrixND<Point> TransMat;

		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				TransMat.m[i1*dim+i2] = m[i2*dim+i1];
			}
		}
		return TransMat;
	}

// ****************************************************************************************

	//! Transposes the matrix in-place.
	inline void transpose() {
		ValueType m2[numFields];
		for (UInt i1 = 0; i1 < numFields; i1++) m2[i1] = m[i1];

		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				m[i1*dim+i2] = m2[i2*dim+i1];
			}
		}

	}

// ****************************************************************************************

	//! Backsolves the SVD-decomposition obtained by 'SVD_decomp(...)'.
	void SVD_bksolve(SquareMatrixND<Point>& U_m, Point& sigma, SquareMatrixND<Point>& V_m, Point& b, Point& x_result) {

		// Build matrices and vectors NR-style:
		ValueType** U_m_NR = new ValueType*[dim+1];
		ValueType** V_m_NR = new ValueType*[dim+1];
		ValueType* sigma_NR = new ValueType[dim+1];
		ValueType* b_NR = new ValueType[dim+1];
		ValueType* x_result_NR = new ValueType[dim+1];

		for(int i1 = 0; i1 < dim+1; i1++) {
			U_m_NR[i1] = new ValueType[dim+1];
			V_m_NR[i1] = new ValueType[dim+1];
		}

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				U_m_NR[i+1][j+1] = U_m(i,j);
				V_m_NR[i+1][j+1] = V_m(i,j);
			}
			sigma_NR[i+1] = sigma[i];
			b_NR[i+1] = b[i];
		}

		svbksb<ValueType>(U_m_NR, sigma_NR, V_m_NR, dim, dim, b_NR, x_result_NR);

		// copy solution:
		for (int i1 = 1; i1 < dim+1; i1++) x_result[i1-1] = x_result_NR[i1];

		// Clean up NR-mess
		for(int i1 = 0; i1 < dim+1; i1++) {
			delete[] U_m_NR[i1];
			delete[] V_m_NR[i1];
		}
		delete[] U_m_NR;
		delete[] V_m_NR;

		delete[] x_result_NR;
		delete[] b_NR;
		delete[] sigma_NR;

	}

// ****************************************************************************************

	//! Computes the SVD-decomposition of the matrix.
	bool SVD_decomp(SquareMatrixND<Point>& U_m, Point& sigma, SquareMatrixND<Point>& V_m) const {

		// Build dim x dim matrix NR-style:
		ValueType** CV = new ValueType*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++) CV[i1] = new ValueType[dim+1];
		for (int i = 0; i < dim; i++) for (int j = 0; j < dim; j++) CV[i+1][j+1] = (*this)(i,j);

        ValueType sigma1[dim+1];
		ValueType** V = new ValueType*[dim+1];
		for(int i1 = 0; i1 < dim+1; i1++) V[i1] = new ValueType[dim+1];

		bool res = svdcmp<ValueType>(CV, dim, dim, sigma1, V);

		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				U_m(i1,i2) = CV[i1+1][i2+1];
				V_m(i1,i2) = V[i1+1][i2+1];
			}
			sigma[i1] = sigma1[i1+1];
		}

		for(int i1 = 0; i1 < dim+1; i1++) delete[] V[i1];
		delete[] V;
		for(int i1 = 0; i1 < dim+1; i1++) delete[] CV[i1];
		delete[] CV;

		return res;
	}

// ****************************************************************************************

	//! The matlab \ operator. i.e., A = B.opBackSlash(C) (A=B\C) finds the Matrix A such that BA=C.
	inline SquareMatrixND<Point> opBackSlash(const SquareMatrixND<Point>& C) const {
//		std::cerr << "inOpBackSlash: Rinv:" << std::endl; (*this).getInverted(0.00000000001).print();
		return (*this).getInverted(0.00000000001)*C;
	}

// ****************************************************************************************

	//! The matlab / operator. i.e., A = B.opSlash(C) (A=B/C) finds the Matrix A such that AC=B.
	inline SquareMatrixND<Point> opSlash(const SquareMatrixND<Point>& C) const {
		return (*this)*C.getInverted(0.00000000001);
	}

// ****************************************************************************************

	//! Returns Matrix^T*Vector.
	inline Point TransposedVecTrans(const Point& pt) const {
		Point ret;

		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				ret[i1] += m[i2*dim+i1]*pt[i2];
			}
		}

		return ret;
	}

// ****************************************************************************************

	//! Returns Matrix*Vector.
	inline Point vecTrans(const Point& pt) const {
		Point ret;

		for (UInt i1 = 0; i1 < dim; i1++) {
			for (UInt i2 = 0; i2 < dim; i2++) {
				ret[i1] += m[i1*dim+i2]*pt[i2];
			}
		}

		return ret;
	}

// ****************************************************************************************

	ValueType getDeterminant();

// ****************************************************************************************

	//! Dataarray.
	ValueType m[numFields];

private:

};

	
template <> inline double SquareMatrixND<vec3d>::getDeterminant() {
	return	(m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7]) -
			(m[0]*m[5]*m[7] + m[1]*m[3]*m[8] + m[2]*m[4]*m[6]);
}

template <> inline float SquareMatrixND<vec3f>::getDeterminant() {
	return	(m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7]) -
			(m[0]*m[5]*m[7] + m[1]*m[3]*m[8] + m[2]*m[4]*m[6]);
}

template <> inline double SquareMatrixND<vec2d>::getDeterminant() {
	return	(m[0]*m[3]) - (m[1]*m[2]);
}

template <> inline float SquareMatrixND<vec2f>::getDeterminant() {
	return	(m[0]*m[3]) - (m[1]*m[2]);
}
// ****************************************************************************************
