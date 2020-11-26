#pragma once


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

#ifndef NO_GSL
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#endif


#ifndef NO_GSL
inline void gslVectorPrint(gsl_vector* v, std::ostream& out) {
	out << "Vector: " << std::endl;
	out << "Dimension: " << v->size << std::endl;

	for (int i = 0; i < (int)v->size; i++) {
		out << gsl_vector_get(v,i) << " ";
	}
}

inline void gslMatrixPrint(gsl_matrix* m, std::ostream& out) {
	out << "Matrix: " << std::endl;
	out << "Dimensions: " << m->size1 << " x " << m->size2 << std::endl;

	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			out << std::setprecision(4) << std::fixed << gsl_matrix_get(m,i,j) << " ";
		}
		out << std::endl;
	}
}
#endif




struct PCA_Basis {


	PCA_Basis(const char* filename, int readNumVecs = 0) {
#ifndef NO_GSL
		pseudo = NULL;
#endif
		loadFromFile(filename, readNumVecs);
	}


	PCA_Basis() {
#ifndef NO_GSL
		pseudo = NULL;
#endif
	}

	~PCA_Basis() {
#ifndef NO_GSL
		if (pseudo != NULL)
			gsl_matrix_free(pseudo);
#endif
	}

	void saveToTextFile(const char* filename) {
		std::ofstream fout(filename);

		for (int i1 = 0; i1 < numVecs; i1++) {
			for (int i2 = 0; i2 < sizeOfVecs; i2++) {
				fout << EVecs[i1][i2];
				if (i2 < sizeOfVecs-1) fout << "\t";
			}
			fout << std::endl;
		}
	}

	void saveToBinaryFile(const char* filename) {
		std::ofstream fout(filename, std::ios::binary);
		saveToBinaryStream(fout);
	}


	void loadFromFile(const char* filename, int readNumVecs = 0) {
		std::ifstream fin(filename, std::ios::binary);

		loadFromStream(fin, readNumVecs);
	}

	void saveToBinaryStream(std::ofstream& fout) const {
		int n = (int)EVecs.size();
		if (n == 0) return;
		int m = (int)EVecs[0].size();
		fout.write((char*)&n, sizeof(int));
		fout.write((char*)&m, sizeof(int));
		fout.write((char*)&EVals[0], sizeof(double)*n);
		for (int i = 0; i < n; i++) {
			fout.write((char*)&EVecs[i][0], sizeof(double)*m);
		}
	}


	void loadFromStream(std::ifstream& fin, int readNumVecs = 0, bool outputStatusMessagesToStderr = false) {
		fin.read((char*)&numVecs, sizeof(int));

		int xx = numVecs;
		if (readNumVecs != 0 && numVecs > readNumVecs) numVecs = readNumVecs;

		fin.read((char*)&sizeOfVecs, sizeof(int));
		EVals.resize(numVecs);
		fin.read((char*)&EVals[0], numVecs*sizeof(double));

		if (xx != numVecs) { // Read remaining vectors
			double* tmp = new double[xx];
			fin.read((char*)tmp, (xx-numVecs)*sizeof(double));
			delete[] tmp;
		}

		EVecs.resize(numVecs);
		for (int i1 = 0; i1 < numVecs; i1++) {
			EVecs[i1].resize(sizeOfVecs);
			fin.read((char*)&EVecs[i1][0], sizeOfVecs*sizeof(double));
		}

		//// Divide eigenvectors by eigenvalue:
		//for (int i1 = 0; i1 < numVecs; i1++) {
		//	double EVal = EVals[i1];
		//	for (int i2 = 0; i2 < sizeOfVecs; i2++) {
		//		EVecs[i1][i2] /= EVal;
		//	}
		//}

		if (outputStatusMessagesToStderr) {
			std::cerr << "Read basis with " << numVecs << " basis vectors, each having " << sizeOfVecs << " entries" << std::endl;
		}

		//std::ofstream etxt("eigenvals.txt");
		//std::cerr << "Eigenvalues" << std::endl;
		//for (int i1 = 0; i1 < numVecs; i1++) {
		//	std::cerr << EVals[i1] << std::endl;
		//	etxt << EVals[i1] << std::endl;
		//	//if (i1 >= 49) exit(1);
		//}
	}


	void checkIsOrtho() const {

		for (int i1 = 0; i1 < numVecs; i1++) {
			for (int i2 = 0; i2 < i1; i2++) {
				std::cerr << "x ";
			}
			for (int i2 = i1; i2 < numVecs; i2++) {
				double res = 0;
				for (int i3 = 0; i3 < sizeOfVecs; i3++) {
					res += EVecs[i1][i3]*EVecs[i2][i3];
				}
				if (fabs(res) < 0.00000001) res = 0;
				if (fabs(res-1) < 0.00000001) res = 1;
				std::cerr << res << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl;
		std::cerr << std::endl;
	}


	void print() const {
		for (int i2 = 0; i2 < numVecs; i2++) {
			//for (int i3 = 0; i3 < sizeOfVecs; i3++) {
			//	std::cerr << EVecs[i2][i3] << std::endl;
			//}
			//std::cerr << std::endl;
			//std::cerr << std::endl;
			std::cerr << "EVal[" << i2 << "] = " << EVals[i2] << std::endl;
		}
		std::cerr << sizeOfVecs << std::endl;
	}


	std::vector<double> projectIntoBasis(const std::vector<double>& v) const {
		if (v.size() != sizeOfVecs) {
			std::cerr << "Incompatible vector size in PCA_Basis::projectIntoBasis(...)" << std::endl;
			std::cerr << "vector is " << v.size() << " basis is " << sizeOfVecs << std::endl;
			exit(1);
		}

		std::vector<double> res;
		res.resize(numVecs);

		for (int i1 = 0; i1 < numVecs; i1++) {
			//res[i1] = multColumn(EVecs[i1], v)*EVals[i1];
			res[i1] = multColumn(EVecs[i1], v);
		}

		return res;
	}


	std::vector<double> projectFromBasis(const std::vector<double>& v) const {
		if (v.size() != numVecs) {
			std::cerr << "Incompatible vector size in PCA_Basis::projectFromBasis(...)" << std::endl;
			exit(1);
		}

		std::vector<double> res;
		res.resize(sizeOfVecs, 0);

		for (int i1 = 0; i1 < sizeOfVecs; i1++) {
			for (int i2 = 0; i2 < numVecs; i2++) {
				//res[i1] += v[i2]*EVecs[i2][i1]/EVals[i2];
				res[i1] += v[i2]*EVecs[i2][i1];
			}
		}

		return res;
	}

	void multVarianceIntoBasis() {
		for (int i1 = 0; i1 < sizeOfVecs; i1++) {
			for (int i2 = 0; i2 < numVecs; i2++) {
				//res[i1] += v[i2]*EVecs[i2][i1]/EVals[i2];
				EVecs[i2][i1] *= sqrt(EVals[i2]);
			}
		}
	}

	
#ifndef NO_GSL
	void doSVD(std::vector<bool> boolflag = std::vector<bool>(0)) {

		const int& m = sizeOfVecs;
		const int& n = numVecs;

		gsl_matrix* mA       = gsl_matrix_alloc(m, n);		// matrix A

		if (boolflag.size() > 0) {

			int cnt = 0;
			std::cerr << "zeroing ";
			for(int col = 0; col < m; col++) {
				if (!boolflag[col/3]) {
					cnt++;
				}
			}
			std::cerr << cnt << " / " << boolflag.size()*3 << " rows" << std::endl;

			for(int row = 0; row < n; row++) {
				for(int col = 0; col < m; col++) {
					double val = EVecs[row][col]/EVals[row];
					if (!boolflag[col/3]) { // boolflag[col/3] /3 because we have a boolflag per vertex (which makes 3 components x,y,z).
						val = 0;
					}
					gsl_matrix_set(mA, col, row, val);
				}
			}
		} else {
			for(int row = 0; row < n; row++) {
				for(int col = 0; col < m; col++) {
					double val = EVecs[row][col]/EVals[row];
					gsl_matrix_set(mA, col, row, val);
				}
			}
		}

		//gslMatrixPrint(mA, std::cerr);

		gsl_matrix* mV       = gsl_matrix_alloc(n, n);		// matrix V from svd A = USV^T
		gsl_vector* vs       = gsl_vector_alloc(n);			// singular values 
		gsl_vector* svd_work = gsl_vector_alloc(n);			// work vector required by gsl svd
		gsl_matrix* mSinv    = gsl_matrix_alloc(n, n);		// pseudo inverse matrix of matrix S from svd
		gsl_matrix* mC       = gsl_matrix_alloc(n, n);	    // temporary matrix will hold V * inv(S);
		pseudo   = gsl_matrix_alloc(n, m);	    // temporary matrix will hold V * inv(S);

		int svd_res = gsl_linalg_SV_decomp(mA, mV, vs, svd_work);

		double eps = 0.0000001;

		// create pseudo inverse of S
		gsl_matrix_set_zero(mSinv);
		for(int i = 0; i < numVecs; i++) {
			double s = gsl_vector_get(vs, i);
			//std::cerr << "Singular value " << i << " " << s << "\n";
			s = s <= eps ? 0.0 : 1.0/s;	
			//s = s / (s*s + tikhonov_alpha * tikhonov_alpha);
			gsl_matrix_set(mSinv, i, i, s);
		}

		// C = V * inv(S) 
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mSinv, 0.0, mC);
		// Pseudo = C * U^T
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0,   mC, mA, 0.0, pseudo);

		//gslMatrixPrint(pseudo, std::cerr);

		gsl_matrix_free(mA);
		gsl_matrix_free(mV);     
		gsl_vector_free(svd_work);
		gsl_matrix_free(mSinv);
		gsl_matrix_free(mC);
	}


	std::vector<double> ProjectIntoPCASpaceUseSVD(const std::vector<double>& v, std::vector<bool> boolflag = std::vector<bool>(0)) {

		if (pseudo == NULL) doSVD(boolflag);

		gsl_vector* vb       = gsl_vector_alloc(sizeOfVecs);			// right hand side
		gsl_vector* vx       = gsl_vector_alloc(numVecs);			// right hand side

		for (int i1 = 0; i1 < sizeOfVecs; i1++) {
			gsl_vector_set(vb, i1, v[i1]);
		}

		// solve the system
		// x = V * b;
		gsl_blas_dgemv(CblasNoTrans, 1.0f, pseudo, vb, 0.0f, vx);

		std::vector<double> x;

		for (int i1 = 0; i1 < numVecs; i1++) {
			x.push_back(gsl_vector_get(vx, i1));
		}

		gsl_vector_free(vb);
		gsl_vector_free(vx);

		return x;
	}
#endif


	void othogonalizeUsingGramSchmidt(std::vector<double> ev0) {
		std::vector< std::vector<double> > E;
		E.push_back(ev0);
		for (int i1 = 1; i1 < numVecs; i1++) {
			std::vector<double> working = EVecs[i1];
			for (unsigned int i2 = 0; i2 < E.size(); i2++) {
				const std::vector<double>& v = E[i2];
				double factor = SkalProd(v, working);
				subtractVectorTimesSkalar(working, v, factor);				
			}
			E.push_back(working);
		}
		EVecs = E;
	}


	std::vector< std::vector<double> > EVecs;
	std::vector<double> EVals;
	
#ifndef NO_GSL
	gsl_matrix* pseudo;
#endif

	int numVecs;
	int sizeOfVecs;


private:


	double SkalProd(const std::vector<double>& v0, const std::vector<double>& v1) {
		double res = 0;
		for (unsigned int i1 = 0; i1 < v0.size(); i1++) res += v0[i1]*v1[i1];
		return res;
	}

	void subtractVectorTimesSkalar(std::vector<double>& working, const std::vector<double>& subtrahend, double factor=1) {
		for (unsigned int i1 = 0; i1 < working.size(); i1++) 
			working[i1] -= subtrahend[i1]*factor;
	}


	inline double multColumn(const std::vector<double>& v0, const std::vector<double>& v1) const {
		double res = 0;
		for (int i1 = 0; i1 < (int)v0.size(); i1++) {
			res += v0[i1]*v1[i1];
		}
		return res;
	}


public:

	void loadFromEVectors(const char* file, int mini, int maxi) {

		numVecs = maxi-mini;

		EVals.reserve(numVecs);
		EVecs.reserve(numVecs);

		int i = 0;
		for (int i1 = mini; i1 < maxi; i1++) {

			std::stringstream fn;
			fn << file << i1 << ".bin";

			if (i1 % 10 == 0)
				std::cerr << "Loading basis (" << (int)(100*i1/maxi) << " %)\r";

			std::ifstream fin(fn.str().c_str(), std::ios::binary);
			double EVal;
			fin.read((char*)&EVal, sizeof(double));

			std::cerr << EVal << std::endl;

			fin.read((char*)&sizeOfVecs, sizeof(int));
			std::vector<double> res;
			res.resize(sizeOfVecs);
			fin.read((char*)&res[0], sizeof(double)*sizeOfVecs);
			EVals.push_back(EVal);
			EVecs.push_back(res);
			//for (int i1 = 0; i1 < 10; i1++) std::cerr << res[i1] << std::endl;
			i++;
		}
		std::cerr << std::endl;

	}
};
