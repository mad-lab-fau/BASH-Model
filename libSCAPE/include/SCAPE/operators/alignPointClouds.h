#pragma once

#include "PointCloud.h"
#include "rigid_transform.h"
#include "SimpleMesh.h"


//! Computes the rigid transformation that centers the pc at 0 and orients it such that the first principal direction points into x direction
template <typename T>
RigidTransform<T> PCA_alignPointCloud(std::vector<T>& pc) {
	
	// Compute cog of both point clouds:
	T cog(0,0,0);
	UInt numVertices = pc.size();
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		cog  += pc[i1];
	}
	cog /= numVertices;

	// Build covariance matrix:
	SquareMatrixND<T> COV;
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		T p = pc[i1]-cog;
		COV.addFromTensorProduct(p, p);
	}

	RigidTransform<T> trans;

	COV.calcEValuesAndVectorsCORRECT(trans.rotate);
	trans.rotate.transpose();
	trans.translate = trans.rotate.vecTrans(-cog);
	return trans;
}

//! Computes and returns the best fitting rigid transformation between the point sets 'pc1' and 'pc2' (i.e., how to rotate/move pc2 such that it matches the alignment of pc1).
template <typename T>
static RigidTransform<T> rigidAlignPointCloudsWeighted(const T* pc1, const T* pc2, UInt numVertices, const typename T::ValueType* weights) {

	typedef T Point; 
	typedef Point::ValueType ValueType;

	// Compute cog of both point clouds:
	Point cog_f_i(0,0,0);
	Point cog_f_i1(0,0,0);
	ValueType sum = 0;
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		cog_f_i  += pc1[i1]*weights[i1];
		cog_f_i1 += pc2[i1]*weights[i1];
		sum += weights[i1];
	}
	cog_f_i /= sum;
	cog_f_i1 /= sum;

	// Build cross-covariance matrix:
	SquareMatrixND<Point> COV;
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		COV.addFromTensorProduct((pc1[i1]-cog_f_i)*weights[i1], (pc2[i1]-cog_f_i1)*weights[i1]);
	}

	// Compute aligning rotation
	SquareMatrixND<Point> U, V;
	Point sigma;
	bool SVD_result = COV.SVD_decomp(U, sigma, V);
	V.transpose();
	SquareMatrixND<Point> rot_matrix = U*V;


	ValueType det = rot_matrix.getDeterminant();
	if (det < 0) {
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
		rot_matrix = U*V;
	}
	//rot_matrix.transpose();

	Point new_cog = rot_matrix.vecTrans(cog_f_i1);
	vec3f cog_offset = cog_f_i - new_cog;

	RigidTransform<Point> trans(cog_offset, rot_matrix);
	return trans;
}


//! Computes and returns the best fitting rigid transformation between the point sets 'pc1' and 'pc2' (i.e., how to rotate/move pc2 such that it matches the alignment of pc1).
template <typename T>
static RigidTransform<T> rigidAlignPointClouds(const T* pc1, const T* pc2, UInt numVertices) {

	typedef T Point; 
	typedef Point::ValueType ValueType;

	// Compute cog of both point clouds:
	Point cog_f_i(0,0,0);
	Point cog_f_i1(0,0,0);
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		cog_f_i  += pc1[i1];
		cog_f_i1 += pc2[i1];
	}
	cog_f_i /= (ValueType)numVertices;
	cog_f_i1 /= (ValueType)numVertices;

	// Build cross-covariance matrix:
	SquareMatrixND<Point> COV;
	for (UInt i1 = 0; i1 < numVertices; i1++) {
		COV.addFromTensorProduct((pc1[i1]-cog_f_i), (pc2[i1]-cog_f_i1));
	}

	// Compute aligning rotation
	SquareMatrixND<Point> U, V;
	Point sigma;
	bool SVD_result = COV.SVD_decomp(U, sigma, V);
	V.transpose();
	SquareMatrixND<Point> rot_matrix = U*V;


	ValueType det = rot_matrix.getDeterminant();
	if (det < 0) {
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
		rot_matrix = U*V;
	}
	//rot_matrix.transpose();

	Point new_cog = rot_matrix.vecTrans(cog_f_i1);
	vec3f cog_offset = cog_f_i - new_cog;

	RigidTransform<Point> trans(cog_offset, rot_matrix);
	return trans;
}

//! Computes and returns the best fitting rigid transformation between the point sets 'pc1' and 'pc2' (i.e., how to rotate/move pc2 such that it matches the alignment of pc1).
template <typename T>
static RigidTransform<T> rigidAlignPointClouds(const T* pc1, const T* pc2, const std::multimap<double, UInt>& indices, UInt numVertices) {

	typedef T Point; 
	typedef Point::ValueType ValueType;

	// Compute cog of both point clouds:
	Point cog_f_i(0,0,0);
	Point cog_f_i1(0,0,0);

	int cnt = 0;
	for (std::multimap<double, UInt>::iterator it = indices.begin(); ((it != indices.end()) && cnt < (numVertices)) ; it++, cnt++) {
		UInt i1 = it->second;
		cog_f_i  += pc1[i1];
		cog_f_i1 += pc2[i1];
	}

	cog_f_i /= (ValueType)numVertices;
	cog_f_i1 /= (ValueType)numVertices;

	// Build cross-covariance matrix:
	SquareMatrixND<Point> COV;
	cnt = 0;
	for (std::multimap<double, UInt>::iterator it = indices.begin(); ((it != indices.end()) && cnt < (numVertices)) ; it++, cnt++) {
		UInt i1 = it->second;
		COV.addFromTensorProduct((pc1[i1]-cog_f_i), (pc2[i1]-cog_f_i1));
	}

	// Compute aligning rotation
	SquareMatrixND<Point> U, V;
	Point sigma;
	bool SVD_result = COV.SVD_decomp(U, sigma, V);
	U.transpose();
	SquareMatrixND<Point> rot_matrix = V*U;
	rot_matrix.transpose();

	Point new_cog = rot_matrix.vecTrans(cog_f_i1);
	vec3f cog_offset = cog_f_i - new_cog;

	RigidTransform<Point> trans(cog_offset, rot_matrix);
	return trans;
}


//! Finds the transformation which aligns 'sm2' to 'sm1'
static RigidTransformWithScale<vec3d> getOptimalAlignment(SimpleMesh* sm1, SimpleMesh* sm2, bool testFlip = true, std::vector<bool> boolflag = std::vector<bool>(0)) {

	typedef vec3d Point; 
	typedef double ValueType;

	UInt numVertices = 0;

	// Compute cog of both point clouds:
	Point cog_f_i(0,0,0);
	Point cog_f_i1(0,0,0);

	if (boolflag.size() == sm1->getNumV()) {
		for (int i1 = 0; i1 < sm1->getNumV(); i1++) {
			if (!boolflag[i1]) continue;
			cog_f_i  += sm1->vList[i1].c;
			cog_f_i1 += sm2->vList[i1].c;
			numVertices++;
		}
	} else {
		numVertices = sm1->getNumV();
		for (UInt i1 = 0; i1 < numVertices; i1++) {
			cog_f_i  += sm1->vList[i1].c;
			cog_f_i1 += sm2->vList[i1].c;
		}
	}

	cog_f_i /= (ValueType)numVertices;
	cog_f_i1 /= (ValueType)numVertices;

	// Build cross-covariance matrix:
	SquareMatrixND<Point> COV;
	double l1 = 0;
	double l2 = 0;

	if (boolflag.size() == sm1->getNumV()) {
		for (int i1 = 0; i1 < sm1->getNumV(); i1++) {
			if (!boolflag[i1]) continue;
			vec3d e1 = (sm1->vList[i1].c-cog_f_i);
			vec3d e2 = (sm2->vList[i1].c-cog_f_i1);
			COV.addFromTensorProduct(e1, e2);
			l1 += e1.length();
			l2 += e2.length();
		}
	} else {
		for (UInt i1 = 0; i1 < numVertices; i1++) {
			vec3d e1 = (sm1->vList[i1].c-cog_f_i);
			vec3d e2 = (sm2->vList[i1].c-cog_f_i1);
			COV.addFromTensorProduct(e1, e2);
			l1 += e1.length();
			l2 += e2.length();
		}
	}

	COV *= (1.0/numVertices);

	// Compute aligning rotation
	SquareMatrixND<Point> U, V;
	Point sigma;
	bool SVD_result = COV.SVD_decomp(U, sigma, V);
	V.transpose();
	SquareMatrixND<Point> rot_matrix = U*V;


	ValueType det = rot_matrix.getDeterminant();
	if ((testFlip) && (det < 0)) {
		std::cerr << "Flip" << std::endl;
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
		rot_matrix = U*V;
	}

	Point new_cog = rot_matrix.vecTrans(cog_f_i1 * l1 / l2);
	vec3f cog_offset = cog_f_i-new_cog;


	//std::cerr << "Scale is: " << l1/l2 << std::endl;

	RigidTransformWithScale<vec3d> rot(cog_offset, rot_matrix, l1/l2);

	return rot;
}


//! Finds the transformation which aligns 'pc2' to 'pc1'
template <class T>
static RigidTransformWithScale<T> getOptimalAlignment(const std::vector<T>& pc1, const std::vector<T>& pc2, bool testFlip = true, std::vector<bool> boolflag = std::vector<bool>(0)) {

	typedef T Point; 
	typedef T::ValueType ValueType;

	UInt numVertices = 0;

	// Compute cog of both point clouds:
	Point cog_f_i(0,0,0);
	Point cog_f_i1(0,0,0);

	if (boolflag.size() == pc1.size()) {
		for (int i1 = 0; i1 < (int)pc1.size(); i1++) {
			if (!boolflag[i1]) continue;
			cog_f_i  += pc1[i1];
			cog_f_i1 += pc2[i1];
			numVertices++;
		}
	} else {
		numVertices = (UInt)pc1.size();
		for (UInt i1 = 0; i1 < numVertices; i1++) {
			cog_f_i  += pc1[i1];
			cog_f_i1 += pc2[i1];
		}
	}

	cog_f_i /= (ValueType)numVertices;
	cog_f_i1 /= (ValueType)numVertices;

	// Build cross-covariance matrix:
	SquareMatrixND<Point> COV;
	ValueType l1 = 0;
	ValueType l2 = 0;

	if (boolflag.size() == numVertices) {
		for (int i1 = 0; i1 < (int)numVertices; i1++) {
			if (!boolflag[i1]) continue;
			vec3d e1 = (pc1[i1]-cog_f_i);
			vec3d e2 = (pc2[i1]-cog_f_i1);
			COV.addFromTensorProduct(e1, e2);
			l1 += (ValueType)e1.length();
			l2 += (ValueType)e2.length();
		}
	} else {
		for (UInt i1 = 0; i1 < numVertices; i1++) {
			vec3d e1 = (pc1[i1]-cog_f_i);
			vec3d e2 = (pc2[i1]-cog_f_i1);
			COV.addFromTensorProduct(e1, e2);
			l1 += (ValueType)e1.length();
			l2 += (ValueType)e2.length();
		}
	}


	//COV *= (1.0/numVertices);

	// Compute aligning rotation
	SquareMatrixND<Point> U, V;
	Point sigma;
	bool SVD_result = COV.SVD_decomp(U, sigma, V);
	V.transpose();
	SquareMatrixND<Point> rot_matrix = U*V;


	ValueType det = rot_matrix.getDeterminant();
	if ((testFlip) && (det < 0)) {
		std::cerr << "Flip" << std::endl;
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
		rot_matrix = U*V;
	}

	Point new_cog = rot_matrix.vecTrans(cog_f_i1*l1/l2);
	Point cog_offset = cog_f_i-new_cog;


	//std::cerr << "Scaler is: " << l1/l2 << std::endl;

	RigidTransformWithScale<Point> rot(cog_offset, rot_matrix, l1/l2);

	return rot;
}
