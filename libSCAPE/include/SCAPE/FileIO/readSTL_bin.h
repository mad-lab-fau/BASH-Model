#pragma once


#include <fstream>
#include <iostream>
#include <vector>

#include "../point3d.h"
#include "../FindDuplicatesOctree.h"

template <class T>
void readBinarySTL(const char* fn, std::vector<T>& pts) {

	std::ifstream myFile (fn, std::ios::in | std::ios::binary);

	char header_info[80] = "";
	unsigned int nTri;
 
	char tmp[2];
	vec3f n, p0, p1, p2;

	//read 80 byte header
	if (myFile) {
		myFile.read(header_info, 80);
		myFile.read((char*)&nTri, 4);
		std::cerr << "numTris: " << nTri << std::endl;

		for (unsigned int i1 = 0; i1 < nTri; i1++) {
			myFile.read((char*)&n, 12);
			myFile.read((char*)&p0, 12);
			myFile.read((char*)&p1, 12);
			myFile.read((char*)&p2, 12);
			myFile.read(tmp, 2);

			pts.push_back(p0);
			pts.push_back(p1);
			pts.push_back(p2);
		}
	} else {
		std::cerr << "Error reading file " << fn << std::endl;
		std::cerr << "in readBinarySTL(char*)" << std::endl;
		exit(1);
	}
}


template <class T>
void readBinarySTLCleanUp(const char* fn, std::vector<T>& verts, std::vector<vec3i>& tris) {

	std::vector<T> pts;
	readBinarySTL(fn, pts);

	T mmin, mmax;
	computeAABB(pts, mmin, mmax);
	T diag = mmax-mmin;
	mmin -= diag*0.01;
	mmax += diag*0.01;


	JBSlib::DFOctree3d octree(10, mmin, mmax);

	UInt pos_p0, pos_p1, pos_p2;

	unsigned int numT = pts.size()/3;
	for (unsigned int i1 = 0; i1 < numT; i1++) {
		const T& p0 = pts[i1*3+0];
		const T& p1 = pts[i1*3+1];
		const T& p2 = pts[i1*3+2];

		if (octree.insertObject(pos_p0, p0)) verts.push_back(p0);
		if (octree.insertObject(pos_p1, p1)) verts.push_back(p1);
		if (octree.insertObject(pos_p2, p2)) verts.push_back(p2);

		tris.push_back(vec3i(pos_p0, pos_p1, pos_p2));
	}

}
