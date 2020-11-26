#pragma once

#include <set>
#include <list>
#include <vector>
#include <iomanip>


// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************


template <class T>
class PairEnt {

public:
	PairEnt(int pos_, T value_) : pos(pos_), value(value_) {}

	bool operator< (const PairEnt& other) const {
		return (pos < other.pos);
	}

	int pos;
	T value;
};


// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************


template <class T>
class RandAccessCompRowMatrix {

public:
	
// ****************************************************************************************

	typedef PairEnt<T> ENT;
	
// ****************************************************************************************

	RandAccessCompRowMatrix(int dim_n, int dim_m) : m_dim_n(dim_n), m_dim_m(dim_m) {
		m_rows.resize(dim_n);
	}
	
// ****************************************************************************************

	void add(int x, int y, T value) {

#ifdef DEBUG_BOUNDS_CHECK
		if ((x >= m_dim_n) || (y >= m_dim_m)) {
			std::cerr << "RandAccessCompRowMatrix: out of bounds access" << std::endl;
			std::cerr << "Accessing " << x << ", " << y << std::endl;
			std::cerr << "MatrixDimensions " << m_dim_n << ", " << m_dim_m << std::endl;
			exit(1);
		}
#endif
		typename std::list<ENT>::iterator it = m_rows[x].begin();

		for (; it != m_rows[x].end() && it->pos < y; it++) {}

		if (it == m_rows[x].end()) {
			ENT e(y, value);
			m_rows[x].insert(it, e);
			return;
		}

		if (it->pos == y) {
			it->value += value;
		} else {
			ENT e(y, value);
			m_rows[x].insert(it, e);
		}
	}

// ****************************************************************************************

	void printRow(int i1, std::ostream& s, int precision = 5) const {
		for (typename std::list<ENT>::const_iterator it = m_rows[i1].begin(); it != m_rows[i1].end(); it++) {
			s << std::fixed << std::setprecision(precision) << std::setw(precision) << it->pos << ": " << it->value << " ";
		}
		s << std::endl;
	}

// ****************************************************************************************

	void print(std::ostream& s, int precision = 5) const {
		for (int i1 = 0; i1 < (int)m_rows.size() && i1 < 10; i1++) {
			int i2 = 0;
			for (typename std::list<ENT>::const_iterator it = m_rows[i1].begin(); it != m_rows[i1].end(); it++) {
				while (i2 < it->pos) {
					i2++;
					s << std::fixed << std::setprecision(precision) << std::setw(precision) << 0 << " ";
				}
				s << std::fixed << std::setprecision(precision) << std::setw(precision) << it->value << " ";
				i2++;
			}
			while (i2 < m_dim_m) {
				i2++;
				s << std::fixed << std::setprecision(precision) << std::setw(precision) << 0 << " ";
			}
			s << std::endl;
		}
	}
	
// ****************************************************************************************

	void getMatrix(std::vector<double>& entries, std::vector<int>& row_index, std::vector<int>& col_ptr) {

		entries.clear();
		row_index.clear();
		col_ptr.clear();
		
		col_ptr.push_back(0);
		for (int x = 0; x < m_dim_n; x++) {

			typename std::list<ENT>::iterator it = m_rows[x].begin();
			for (; it != m_rows[x].end(); it++) {
				entries.push_back(it->value);
				row_index.push_back(it->pos);
			}

			col_ptr.push_back(entries.size());
		}
	}

// ****************************************************************************************

private:

	std::vector<std::list<ENT>> m_rows;

	//! height of the matrix.
	int m_dim_n;

	//! width of the matrix.
	int m_dim_m;

};




// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************





template <class T>
class RandAccessCompColMatrix {

public:
	
// ****************************************************************************************

	typedef PairEnt<T> ENT;

// ****************************************************************************************

	RandAccessCompColMatrix(int dim_n, int dim_m) : m_dim_n(dim_n), m_dim_m(dim_m) {
		m_cols.resize(dim_m);
	}

// ****************************************************************************************

	void add(int x, int y, T value) {

#ifdef DEBUG_BOUNDS_CHECK
		if ((x >= m_dim_n) || (y >= m_dim_m)) {
			std::cerr << "RandAccessCompRowMatrix: out of bounds access" << std::endl;
			std::cerr << "Accessing " << x << ", " << y << std::endl;
			std::cerr << "MatrixDimensions " << m_dim_n << ", " << m_dim_m << std::endl;
			exit(1);
		}
#endif

		typename std::list<ENT>::iterator it = m_cols[y].begin();

		for (; it != m_cols[y].end() && it->pos < x; it++) {}

		if (it == m_cols[y].end()) {
			ENT e(x, value);
			m_cols[y].insert(it, e);
			return;
		}

		if (it->pos == x) {
			it->value += value;
		} else {
			ENT e(x, value);
			m_cols[y].insert(it, e);
		}
	}

// ****************************************************************************************

	std::vector<T> getAtb(std::vector<T>& b) {

#ifdef DEBUG_BOUNDS_CHECK
		if ((int)b.size() != m_dim_n) {
			std::cerr << "RandAccessCompRowMatrix: incompatible vector size in getAtb()" << std::endl;
			exit(1);
		}
#endif

		std::vector<T> Atb;
		Atb.reserve(m_dim_m);

		for (int col = 0; col < m_dim_m; col++) {
			T value = 0;
			for (typename std::list<ENT>::iterator it = m_cols[col].begin(); it != m_cols[col].end(); it++) {
				value += it->value*b[it->pos];
			}
			Atb.push_back(value);
		}

		return Atb;

	}

// ****************************************************************************************

	void setBottomDown(int x, int y, T value) {

#ifdef DEBUG_BOUNDS_CHECK
		if ((x >= m_dim_n) || (y >= m_dim_m)) {
			std::cerr << "RandAccessCompRowMatrix: out of bounds access" << std::endl;
			std::cerr << "Accessing " << x << ", " << y << std::endl;
			std::cerr << "MatrixDimensions " << m_dim_n << ", " << m_dim_m << std::endl;
			exit(1);
		}
#endif

		//! don't save zeros.
		if (value == 0) return;

		ENT e(x, value);
		m_cols[y].push_back(e);
	}

// ****************************************************************************************

	void getAtAForUMFPACKPrecomputeNonZeroEntries(std::vector<T>& entries, std::vector<int>& row_index, std::vector<int>& col_ptr) {

		// Convert to compressed row matrix
		std::vector<std::vector<int>> rows;
		rows.resize(m_dim_n);
		for (int col = 0; col < m_dim_m; col++) {
			for (typename std::list<ENT>::iterator it = m_cols[col].begin(); it != m_cols[col].end(); it++) {
				rows[it->pos].push_back(col);
			}
		}

		std::vector<std::set<int>> neighbors;
		neighbors.resize(m_dim_m);
		for (int row = 0; row < m_dim_n; row++) {
			for (int col = 0; col < (int)rows[row].size(); col++) {
				const int& v0 = rows[row][col];
				for (int col2 = 0; col2 < (int)rows[row].size(); col2++) {
					const int& v1 = rows[row][col2];
					neighbors[v0].insert(v1);
				}
			}
		}


		col_ptr.push_back(0);
		for (int x = 0; x < m_dim_m; x++) {

			for (std::set<int>::iterator it = neighbors[x].begin(); it != neighbors[x].end(); ++it) {
				const int& n = *it;
				T val = multRows(x,n);
				entries.push_back(val);
				row_index.push_back(n);
			}

			col_ptr.push_back((int)entries.size());
		}

	}	

// ****************************************************************************************

	void getAtAForTAUCSPrecomputeNonZeroEntries(std::vector<T>& entries, std::vector<int>& row_index, std::vector<int>& col_ptr) {

		// Convert to compressed row matrix
		std::vector<std::vector<int>> rows;
		rows.resize(m_dim_n);
		for (int col = 0; col < m_dim_m; col++) {
			for (typename std::list<ENT>::iterator it = m_cols[col].begin(); it != m_cols[col].end(); it++) {
				rows[it->pos].push_back(col);
			}
		}

		std::vector<std::set<int>> neighbors;
		neighbors.resize(m_dim_m);
		for (int row = 0; row < m_dim_n; row++) {
			for (int col = 0; col < (int)rows[row].size(); col++) {
				const int& v0 = rows[row][col];
				for (int col2 = 0; col2 < (int)rows[row].size(); col2++) {
					const int& v1 = rows[row][col2];
					if (v1 < v0) continue;
					neighbors[v0].insert(v1);
				}
			}
		}


		col_ptr.push_back(0);
		for (int x = 0; x < m_dim_m; x++) {

			for (std::set<int>::iterator it = neighbors[x].begin(); it != neighbors[x].end(); ++it) {
				const int& n = *it;
				T val = multRows(x,n);
				entries.push_back(val);
				row_index.push_back(n);
			}

			col_ptr.push_back((int)entries.size());
		}

	}	

// ****************************************************************************************

	void getAtAForTAUCS(std::vector<T>& entries, std::vector<int>& row_index, std::vector<int>& col_ptr) {

		col_ptr.push_back(0);

		for (int x = 0; x < m_dim_m; x++) {

			for (int y = x; y < m_dim_m; y++) {
				T val = multRows(x,y);
				if (val != 0) {
					//std::cerr << "(" << x+1 << "," << y+1 << ") \t " << val << std::endl;
					entries.push_back(val);
					row_index.push_back(y);
				}
			}

			col_ptr.push_back(entries.size());
		}

	}
	
// ****************************************************************************************

	void print() {
		std::cerr << "Rows: " << m_dim_n << " \t Cols: " << m_dim_m << std::endl;
		std::cerr << "Caution, representation starts at (1,1) for MatLab compatibility!" << std::endl;
		for (int col = 0; col < m_dim_m; col++) {
			for (typename std::list<ENT>::iterator it = m_cols[col].begin(); it != m_cols[col].end(); it++) {
				std::cerr << "(" << it->pos+1 << "," << col+1 << ") \t " << it->value << std::endl;
			}
		}
	}

// ****************************************************************************************

private:

	T multRows(int i, int j) {
		T val = 0;

		typename std::list<ENT>::iterator it_i = m_cols[i].begin();
		typename std::list<ENT>::iterator it_j = m_cols[j].begin();

		while (it_i != m_cols[i].end() && it_j != m_cols[j].end()) {

			//std::cerr << it_i->pos << " " << it_j->pos << " \t " << val << std::endl;

			if (it_i->pos < it_j->pos) it_i++;
			else if (it_i->pos > it_j->pos) it_j++;
			else {
				val += it_i->value*it_j->value;
				it_i++;
				it_j++;
			}
		}

		return val;
	}


	std::vector<std::list<ENT>> m_cols;

	//! height of the matrix.
	int m_dim_n;

	//! width of the matrix.
	int m_dim_m;

};

