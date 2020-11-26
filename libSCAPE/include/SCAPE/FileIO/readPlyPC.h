#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

#include "../SimpleMesh.h"
#include "../point3d.h"

struct plyFile {

	UInt num_cols;
	UInt num_rows;
	UInt numV;
	UInt numG;
	UInt numF;

	plyFile() : num_cols(0), num_rows(0), numV(0), numG(0), numF(0), vertices(0), tris(0), pixmap(0) {}

	plyFile(const plyFile& o) : num_cols(o.num_cols), num_rows(o.num_rows), numV(o.numV), numG(o.numG), numF(o.numF), vertices(0), tris(0), pixmap(0) {

		if (numV > 0) {
			vertices = new vec3d[numV];
			for (UInt i1 = 0; i1 < numV; i1++) vertices[i1] = o.vertices[i1];
		}
		if (numG > 0) {
			pixmap = new UInt[numG];
			for (UInt i1 = 0; i1 < numG; i1++) pixmap[i1] = o.pixmap[i1];
		}
		if (numF > 0) {
			tris = new vec3ui[numF];
			for (UInt i1 = 0; i1 < numF; i1++) tris[i1] = o.tris[i1];
		}

	}

	~plyFile() {
		if (vertices != NULL) {
			delete[] vertices;
			vertices = NULL;
		}
		if (tris != NULL) {
			delete[] tris;
			tris = NULL;
		}
		if (pixmap != NULL) {
			delete[] pixmap;
			pixmap = NULL;
		}
	}

	void setHeader(std::string s) {
		char first_word[100];

		std::stringstream ss(s);
		ss >> first_word;
		std::string fw(first_word);

		std::cerr << s << std::endl;


		if (fw.compare(std::string("obj_info")) == 0) {
			char name[100];
			ss >> name;
			std::string n(name);

			//std::cerr << "obj_info: " << n << std::endl;

			int ii;
			//double dd;

			if (n.compare(std::string("num_cols")) == 0) {
				ss >> ii;
				num_cols = ii;
			} else if (n.compare(std::string("num_rows")) == 0) {
				ss >> ii;
				num_rows = ii;
			}
		}
		if (fw.compare(std::string("element")) == 0) {
			char name[100];
			ss >> name;
			std::string n(name);

			//std::cerr << "element: " << n << std::endl;

			int ii;
			//double dd;

			if (n.compare(std::string("vertex")) == 0) {
				ss >> ii;
				numV = ii;
			} else if (n.compare(std::string("range_grid")) == 0) {
				ss >> ii;
				numG = ii;
			} else if (n.compare(std::string("face")) == 0) {
				ss >> ii;
				numF = ii;
			}
		}
	}

	vec3d* vertices;
	vec3ui* tris;
	UInt* pixmap;
};

plyFile parsePlyFileAscii(const char* filename) {

	plyFile ret;

	std::ifstream file(filename);

	if (!file.good()) return ret;

	std::string s;

	// read header

	std::getline(file, s);
	int i1 = 0;
	while (s.compare(std::string("end_header"))) {
		ret.setHeader(s);
		std::getline(file, s);
	}

	ret.vertices = new vec3d[ret.numV];
	double x,y,z,tmp;
	tmp = 0;
	for (UInt i1 = 0; i1 < ret.numV; i1++) {
		//file >> x >> y >> z >> tmp >> tmp >> tmp;
		file >> x >> y >> z;
		vec3d pt(x,y,z);
		ret.vertices[i1] = pt;
	}


	if (ret.numG > 0) {
		UInt i;
		ret.pixmap = new UInt[ret.numG];
		for (UInt i1 = 0; i1 < ret.numG; i1++) {
			file >> i;

			int ind;

			if (i == 1) {
				file >> ind;
				ret.pixmap[i1] = ind;
			} else {
				ret.pixmap[i1] = 0;
			}

		}

		std::cerr << " numV: " << ret.numV << " numG: " << ret.numG << std::endl;

	} else if (ret.numF > 0) {

		std::cerr << "reading triangles" << std::endl;

		UInt i, v0, v1, v2;
		ret.tris = new vec3ui[ret.numF];
		for (UInt i1 = 0; i1 < ret.numF; i1++) {
			file >> i >> v0 >> v1 >> v2;
			//std::cerr << i << " " << v0 << " " << v1 << " " << v2 << std::endl;
			ret.tris[i1] = vec3ui(v0, v1, v2);
		}

	}

	return ret;
}

SimpleMesh* readPlyFileAscii(const char* filename, double depthdisp = 0) {

	plyFile ply = parsePlyFileAscii(filename);


	SimpleMesh* sm = new SimpleMesh;
	UInt dx = ply.num_cols;
	UInt dy = ply.num_rows;
	for (UInt i = 0; i < ply.numV; i++) {
		sm->insertVertex(ply.vertices[i]);
	}

	if (ply.numG > 0) {

		int i1, i2, i3, i4;

		if (depthdisp == 0) depthdisp = std::numeric_limits<double>::max();

		for (UInt x = 0; x < (dx-1); x++) {
			for (UInt y = 0; y < (dy-1); y++) {
				i1 = ply.pixmap[(x+0)+(y+0)*dx];
				i2 = ply.pixmap[(x+1)+(y+0)*dx];
				i3 = ply.pixmap[(x+0)+(y+1)*dx];
				i4 = ply.pixmap[(x+1)+(y+1)*dx];

				double depthdisp1 = std::max(std::max(fabs(ply.vertices[i1].z-ply.vertices[i3].z), fabs(ply.vertices[i1].z-ply.vertices[i4].z)), fabs(ply.vertices[i3].z-ply.vertices[i4].z));
				double depthdisp2 = std::max(std::max(fabs(ply.vertices[i1].z-ply.vertices[i2].z), fabs(ply.vertices[i1].z-ply.vertices[i4].z)), fabs(ply.vertices[i2].z-ply.vertices[i4].z));

				if (i1 > 0 && i2 > 0 && i3 > 0 && i4 > 0) {
					if (depthdisp1 < depthdisp) sm->insertTriangle(i1, i4, i3);
					if (depthdisp2 < depthdisp) sm->insertTriangle(i4, i1, i2);
				}
			}
		}
	}

	if (ply.numF > 0) {

		for (UInt i1 = 0; i1 < ply.numF; i1++) {
			const vec3ui& t = ply.tris[i1];
			sm->insertTriangle(t.x, t.y, t.z);
		}

	}

	return sm;
}


template <class T>
void readPlyFileAscii(std::vector< T >& verts, std::vector< vec3i >& tris, const char* filename) {

	plyFile ply = parsePlyFileAscii(filename);


	for (UInt i = 0; i < ply.numV; i++) {
		verts.push_back(ply.vertices[i]);
	}

	if (ply.numF > 0) {
		for (UInt i1 = 0; i1 < ply.numF; i1++) {
			const vec3ui& t = ply.tris[i1];
			tris.push_back(t);
		}
	}

}
