#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "point3d.h"
#include "rigid_transform.h"

typedef unsigned int UInt;

struct Part {

	typedef unsigned int UInt;

	int partID;
	std::vector<UInt> verts;
	std::vector<UInt> tris;
	RigidTransform<vec3d> R;

	void savePart(std::ofstream& fout) const {
		fout.write((char*)&partID, sizeof(int));
		int numV = (int)verts.size();
		fout.write((char*)&numV, sizeof(int));
		fout.write((char*)&verts[0], sizeof(UInt)*numV);
		int numT = (int)tris.size();
		fout.write((char*)&numT, sizeof(int));
		fout.write((char*)&tris[0], sizeof(UInt)*numT);
		fout.write((char*)&R, sizeof(RigidTransform<vec3d>));
	}

	void readPart(std::ifstream& fin) {
		fin.read((char*)&partID, sizeof(int));
		int numV;
		fin.read((char*)&numV, sizeof(int));
		verts.clear();
		verts.resize(numV);
		fin.read((char*)&verts[0], sizeof(UInt)*numV);
		int numT;
		fin.read((char*)&numT, sizeof(int));
		tris.clear();
		tris.resize(numT);
		fin.read((char*)&tris[0], sizeof(UInt)*numT);
		fin.read((char*)&R, sizeof(RigidTransform<vec3d>));
	}
};


struct Joint {
	int jointID;
	int part0;
	int part1;
	vec3d pos;
	vec3d rot;

	void print() const {
		std::cerr << jointID << " [" << part0 << " - " << part1 << " \t";
		pos.print();
	}
};


static std::vector<Part> readParts(const char* fn, std::vector<UInt>& lookup) {

	std::ifstream fin(fn);
	UInt numV;
	int numParts;
	fin >> numV >> numParts;

	std::vector<Part> parts(numParts);
	std::cerr << parts.size() << std::endl;

	int p;
	for (UInt i1 = 0; i1 < numV; i1++) {
		fin >> p;
		lookup.push_back(p);
		parts[p].verts.push_back(i1);
	}


	return parts;
}


static std::vector<Joint> readJoints(const char* fn, std::vector< std::vector <int> >& partToJoint) {

	std::vector<Joint> joints;

	std::ifstream fin(fn);
	int numJoints;
	fin >> numJoints;

	partToJoint.resize(20);

	for (int i1 = 0; i1 < numJoints; i1++) {
		Joint j;
		j.jointID = i1;
		fin >> j.part0 >> j.part1 >> j.pos.x >> j.pos.y >> j.pos.z;
		joints.push_back(j);
		partToJoint[j.part0].push_back(i1);
		partToJoint[j.part1].push_back(i1);
	}

	return joints;
}



struct perTriangleParams {

	double vals[3][3][7];

};


struct EdgePair {
	vec3d v0, v1;
	EdgePair() : v0(0.0), v1(0.0) {}
	EdgePair(const vec3d& v0_, const vec3d& v1_) : v0(v0_), v1(v1_) {}
};
