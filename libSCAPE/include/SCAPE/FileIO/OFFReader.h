#pragma once


#include <vector>
#include <fstream>
#include <string>

#include "../point3d.h"

//! Reads the vertices and the triangles of the OFF file 'filename' into the vectors 'verts' and 'tris'. Returns 'true' on success.
template <class T>
bool readOFFFile(std::vector< T >& verts, std::vector< point3d<int> >& tris, const char* filename) {

	verts.clear();
	tris.clear();

	std::ifstream InFile(filename);
	// check if open succeeded
	if (InFile.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	// first line should say 'OFF'
	char string1[5];
	InFile >> string1;

	// read header
	int numV = 0;
	int numT = 0;
	int numE = 0;
	InFile >> numV >> numT >> numE;

	std::cerr << "reading file: " << filename << " " << numV << " verts, " << numT << " tris" << std::endl;

	verts.reserve(numV);
	tris.reserve(numT);

	// read points
	T pt;
	for(int i1 = 0; i1 < numV; i1++) {
		InFile >> pt.x >> pt.y >> pt.z;
		verts.push_back(pt);
	}

	// read triangles
	int num_vs;
	vec3i t;
	for(int i1 = 0; i1 < numT; i1++) {
		InFile >> num_vs;
		if (num_vs > 4) {
			std::cerr << "error reading file " << filename << std::endl << "program only supports triangular and quad meshes" << std::endl;
			return false;
		}
		InFile >> t.x >> t.y >> t.z;
		tris.push_back(t);
		if (num_vs == 4) {
			vec3i t2(t.x, t.z, 0);
			InFile >> t2.z;
			tris.push_back(t2);
		}
	}

	return true;
}


//! Reads the vertices and the triangles and the colors of the COFF file 'filename' into the vectors 'verts' and 'tris' and 'cols'. Returns 'true' on success.
template <class T>
bool readCOFFFile(std::vector< T >& verts, std::vector< point3d<int> >& tris, std::vector<vec3f>& cols, const char* filename) {

	verts.clear();
	tris.clear();
	cols.clear();

	std::ifstream InFile(filename);
	// check if open succeeded
	if (InFile.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	// first line should say 'COFF'
	char string1[5];
	InFile >> string1;

	// read header
	int numV = 0;
	int numT = 0;
	int numE = 0;
	InFile >> numV >> numT >> numE;

	verts.reserve(numV);
	tris.reserve(numT);

	// read points
	T pt;
	vec3f col;
	float tmp;
	for(int i1 = 0; i1 < numV; i1++) {
		InFile >> pt.x >> pt.y >> pt.z;
		InFile >> col.x >> col.y >> col.z >> tmp;
		cols.push_back(col);
		verts.push_back(pt);
	}

	// read triangles
	int num_vs;
	vec3i t;
	for(int i1 = 0; i1 < numT; i1++) {
		InFile >> num_vs;
		if (num_vs != 3) {
			std::cerr << "error reading file " << filename << std::endl << "program only supports triangular meshes" << std::endl;
			return false;
		}
		InFile >> t.x >> t.y >> t.z;
		tris.push_back(t);
	}

	return true;
}


template <class T, class U>
bool writeOFFFile(const std::vector< T >& verts, const std::vector< point3d< U > >& tris, const char* filename) {

	std::ofstream fout(filename);
	// check if open succeeded
	if (fout.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	fout << "OFF\n";
	fout << verts.size() << " " << tris.size() << " 0\n";

	for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
		fout << verts[i1].x << " " << verts[i1].y << " " << verts[i1].z << "\n";
	}
	for (unsigned int i1 = 0; i1 < tris.size(); i1++) {
		fout << "3 " << tris[i1].x << " " << tris[i1].y << " " << tris[i1].z << "\n";
	}

	fout.close();

	return true;

}


template <class T, class U>
bool writeOFFFileWithMarker(const std::vector< T >& verts, const std::vector< point3d< U > >& tris, const char* filename, const std::vector<int>& marker, const std::vector<vec3f>& cols) {

	std::ofstream fout(filename);
	// check if open succeeded
	if (fout.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	if (cols.size() == 0) fout << "OFF\n";
	else fout << "COFF\n";
	fout << verts.size() << " " << tris.size() << " 0\n";

	if (cols.size() == 0) {
		for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
			fout << verts[i1].x << " " << verts[i1].y << " " << verts[i1].z << "\n";
		}
	} else {
		for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
			fout << verts[i1].x << " " << verts[i1].y << " " << verts[i1].z << " "  << cols[i1].x << " " << cols[i1].y << " " << cols[i1].z << " 255\n";
		}
	}

	for (unsigned int i1 = 0; i1 < tris.size(); i1++) {
		fout << "3 " << tris[i1].x << " " << tris[i1].y << " " << tris[i1].z << "\n";
	}

	fout << marker.size() << std::endl;
	for (unsigned int i1 = 0; i1 < marker.size(); i1++) fout << marker[i1] << std::endl;

	fout.close();

	return true;

}


//! Reads the vertices and the triangles of the OFF file 'filename' into the vectors 'verts' and 'tris'. Returns 'true' on success.
template <class T>
bool readOFFFileWithMarker(std::vector< T >& verts, std::vector< point3d<int> >& tris, const char* filename, std::vector<int>& marker, std::vector<vec3f>& colors) {

	verts.clear();
	tris.clear();
	colors.clear();

	std::ifstream InFile(filename);
	// check if open succeeded
	if (InFile.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	// first line should say 'OFF'
	std::string string1;
	InFile >> string1;



	// read header
	int numV = 0;
	int numT = 0;
	int numE = 0;
	InFile >> numV >> numT >> numE;

	verts.reserve(numV);
	tris.reserve(numT);

	// read points
	T pt;


	if (string1.compare("COFF") == 0) {
		colors.reserve(numV);
		vec3f col;
		int t;
		for(int i1 = 0; i1 < numV; i1++) {
			InFile >> pt.x >> pt.y >> pt.z >> col.x >> col.y >> col.z >> t;
			verts.push_back(pt);
			colors.push_back(col);
		}
	} else {
		for(int i1 = 0; i1 < numV; i1++) {
			InFile >> pt.x >> pt.y >> pt.z;
			verts.push_back(pt);
		}
	}

	// read triangles
	int num_vs;
	vec3i t;
	for(int i1 = 0; i1 < numT; i1++) {
		InFile >> num_vs;
		if (num_vs != 3) {
			std::cerr << "error reading file " << filename << std::endl << "program only supports triangular meshes" << std::endl;
			return false;
		}
		InFile >> t.x >> t.y >> t.z;
		tris.push_back(t);
	}

	if (InFile.eof() || !InFile.good()) return true;
	int numMarker;
	InFile >> numMarker;
	if (InFile.eof() || !InFile.good()) return true;
	for (int i1 = 0; i1 < numMarker; i1++) {
		int tmp;
		InFile >> tmp;
		marker.push_back(tmp);
	}

	return true;
}

template <class T>
bool writeCOFFFile(std::vector< T >& verts, std::vector< point3d<int> >& tris, std::vector< point3d<float> >& cols, const char* filename, int multiplier=255) {

	std::ofstream fout(filename);
	// check if open succeeded
	if (fout.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	fout << "COFF\n";
	fout << verts.size() << " " << tris.size() << " 0\n";

	for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
		fout << verts[i1].x << " " << verts[i1].y << " " << verts[i1].z << " " << cols[i1].x*multiplier << " " << cols[i1].y*multiplier << " " << cols[i1].z*multiplier << " " << multiplier << "\n";
	}
	for (unsigned int i1 = 0; i1 < tris.size(); i1++) {
		fout << "3 " << tris[i1].x << " " << tris[i1].y << " " << tris[i1].z << "\n";
	}

	fout.close();

	return true;

}


template <class T>
bool writeCOFFFileNoMult(std::vector< T >& verts, std::vector< point3d<int> >& tris, std::vector< point3d<float> >& cols, const char* filename) {

	std::ofstream fout(filename);
	// check if open succeeded
	if (fout.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	fout << "COFF\n";
	fout << verts.size() << " " << tris.size() << " 0\n";

	for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
		fout << verts[i1].x << " " << verts[i1].y << " " << verts[i1].z << " " << (int)cols[i1].x << " " << (int)cols[i1].y << " " << (int)cols[i1].z << " 255\n";
	}
	for (unsigned int i1 = 0; i1 < tris.size(); i1++) {
		fout << "3 " << tris[i1].x << " " << tris[i1].y << " " << tris[i1].z << "\n";
	}

	fout.close();

	return true;

}

template <class T>
bool writePointCloudOBJ(std::vector< T >& verts, const char* filename) {

	std::ofstream fout(filename);
	// check if open succeeded
	if (fout.fail()) {
		std::cerr << "could not open " << filename << std::endl;
		return false;
	}

	for (unsigned int i1 = 0; i1 < verts.size(); i1++) {
		const T& p = verts[i1];
		fout << "v " << p.x << " " << p.y << " " << p.z << std::endl;
	}
}
