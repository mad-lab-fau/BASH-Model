#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>       // std::toupper
#include <limits>

#include "../SimpleMesh.h"

//! Creates a SimpleMesh from a file.
class MeshReader
{
public:

	//! Constructor.
	MeshReader(){};

	//******************************************************************************************

	static SimpleMesh* readFileFromExtension(const char* filename) {

		SimpleMesh* m_mesh = new SimpleMesh();

		std::string s(filename);
		UInt pos = (UInt) s.rfind(".");
		UInt length = (UInt) s.length();
        
		//std::cerr << "Reading " << filename << std::endl;

		if (pos == s.npos) {
			std::cerr << "ERROR in MeshReaderNew::readFileFromExtension(...)" << std::endl;
			std::cerr << "could not determine file extension from " << s << std::endl;
			std::cerr << "EXIT" << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string ext = s.substr(pos+1);

		std::transform(ext.begin(), ext.end(), ext.begin(), (int(*)(int)) std::tolower);

		if (ext.compare(std::string("off")) == 0) {
			//std::cerr << "Reading off file:" << std::endl;
			readOFFFile(filename, m_mesh);
		} else if (ext.compare(std::string("obj")) == 0) {
			std::cerr << "Reading obj file:" << std::endl;
			readOBJFile(filename, m_mesh);
		} else if (ext.compare(std::string("mesh")) == 0) {
			delete m_mesh;
			std::cerr << "Reading mesh file:" << std::endl;
			m_mesh = readMESHFile(filename);
		} else if (ext.compare(std::string("gts")) == 0) {
			delete m_mesh;
			std::cerr << "Reading gts file:" << std::endl;
			m_mesh = readGTSFile(filename);
		} else if (ext.compare(std::string("bin")) == 0) {
			delete m_mesh;
			std::cerr << "Reading bin file:" << std::endl;
			m_mesh = readBINFile(filename);
		} else if (ext.compare(std::string("cbj")) == 0) {
			delete m_mesh;
			std::cerr << "Reading cbj file:" << std::endl;
			m_mesh = readOBJFile(filename);
		} else if (ext.compare(std::string("cof")) == 0) {
			m_mesh = readColoredOFFFile(filename);
		} else {
			std::cerr << "Unknown extension: " << ext << " specified" << std::endl;
			delete m_mesh;
			m_mesh = NULL;
		}

		return m_mesh;
	}

	//******************************************************************************************

	static SimpleMesh* readBQNFile(const char* filename) {

		std::cerr << "Read quirins binary format" << std::endl;

		std::ifstream fin(filename, std::ifstream::binary);
		if (!fin.good()) {
			std::cerr << "Error reading file " << filename << std::endl;
			return NULL;
		}

		int numV;
		int numT;
		int numA;
		fin.read((char*)&numV, sizeof(int));
		fin.read((char*)&numT, sizeof(int));
		fin.read((char*)&numA, sizeof(int));

		std::cerr << "We have " << numV << " vertices, " << numT << " triangles and " << numA << " attributes" << std::endl;

		std::vector<int> m_attributes;
		m_attributes.resize(numA);
		if(numA!=0){
			fin.read((char*) &m_attributes[0], numA*sizeof(int));	
		}


		vec3f* verts = new vec3f[numV];
		vec3i* tris = new vec3i[numT];

		fin.read((char*)verts, 3*numV*sizeof(float));
		fin.read((char*)tris, 3*numT*sizeof(int));

		SimpleMesh* sm = new SimpleMesh;

		for (int i1 = 0; i1 < numV; i1++) sm->insertVertex(verts[i1].x, verts[i1].y, verts[i1].z);
		for (int i1 = 0; i1 < numT; i1++) sm->insertTriangle(tris[i1].x, tris[i1].y, tris[i1].z);

		delete[] verts;
		delete[] tris;

		std::cerr << "read" << std::endl;

		return sm;
	}

	//******************************************************************************************

	static SimpleMesh* readBINFile(const char* filename) {

		std::ifstream fin(filename, std::ifstream::binary);
		if (!fin.good()) {
			std::cerr << "Error reading file " << filename << std::endl;
			return NULL;
		}

		int numV;
		int numT;
		fin.read((char*)&numV, sizeof(int));
		fin.read((char*)&numT, sizeof(int));

		vec3d* verts = new vec3d[numV];
		vec3i* tris = new vec3i[numT];

		fin.read((char*)verts, 3*numV*sizeof(double));
		fin.read((char*)tris, 3*numT*sizeof(int));


		SimpleMesh* sm = new SimpleMesh;

		for (int i1 = 0; i1 < numV; i1++) sm->insertVertex(verts[i1]);
		for (int i1 = 0; i1 < numT; i1++) sm->insertTriangle(tris[i1].x, tris[i1].y, tris[i1].z);

		delete[] verts;
		delete[] tris;

		return sm;
	}

	//******************************************************************************************

	static void readOFFFile(const char* filename, SimpleMesh* m_mesh) {
		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei!" << std::endl;
			std::cerr << "could not open " << filename << std::endl;
			return;
		}

		char string1[5];
		InFile >> string1;

		int numV = 0;
		int numT = 0;
		int numE = 0;

		InFile >> numV >> numT >> numE;

		// Read points
		double x, y, z;
		for(int i1 = 0; i1 < numV; i1++) {
			InFile >> x >> y >> z;

			//std::cerr << i1 << " " << x << " " << y << " " << z << std::endl;
			m_mesh->insertVertex(x,y,z);
		}

		// Read triangles
		int num_vs, v0, v1, v2, v3, v4, v5;
		for(int i1 = 0; i1 < numT; i1++) {
			InFile >> num_vs;

			//std::cerr << num_vs << std::endl;

			//if (i1%1000 == 0) std::cerr << i1 << std::endl;
			
			if (num_vs == 3) { // A triangle follows:
				InFile >> v0 >> v1 >> v2;
				//std::cerr << "Insert tri " << i1 << ": " << v0 << " " << v1 << " " << v2 << std::endl;
				m_mesh->insertTriangle(v0, v1, v2);
			} else if (num_vs == 4) {// A quadrilataral follows:
				InFile >> v0 >> v1 >> v2 >> v3;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v0, v2, v3);
			} else if (num_vs == 5) {
				InFile >> v0 >> v1 >> v2 >> v3 >> v4;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v2, v3, v4);
				m_mesh->insertTriangle(v2, v4, v0);
			} else if (num_vs == 6) {
				InFile >> v0 >> v1 >> v2 >> v3 >> v4 >> v5;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v2, v3, v4);
				m_mesh->insertTriangle(v2, v4, v5);
				m_mesh->insertTriangle(v2, v5, v0);
			} else if (num_vs > 6) {
				int* vs = new int[num_vs];
				for (int i1 = 0; i1 < num_vs; i1++)
					InFile >> vs[i1];

				m_mesh->insertTriangle(vs[0], vs[1], vs[2]);
				for (int i1 = 0; i1 < num_vs-4; i1++)
					m_mesh->insertTriangle(vs[2], vs[i1+3], vs[i1+4]);

				m_mesh->insertTriangle(vs[2], vs[num_vs-1], vs[0]);
				delete[] vs;
			}
		}

		//std::cerr << "Read " << numV << " Verts " << numT << " Tris" << std::endl;

		//if (numE != 0)
		//	std::cerr << "Warnung bei bool MeshReader::readOFFFile(char *filename)! Kann keine Kanten einlesen!" << std::endl;

		InFile.close();
	}

	//******************************************************************************************

	//! Creates a SimpleMesh from an .OFF file.
	static SimpleMesh* readColoredOFFFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();

		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei!" << std::endl;
			std::cerr << "could not open " << filename << std::endl;
			return m_mesh;
		}

		char string1[5];
		InFile >> string1;

		int numV = 0;
		int numT = 0;
		int numE = 0;

		InFile >> numV >> numT >> numE;

		// Read points
		double x, y, z;
		vec3f col;
		double a;
		for(int i1 = 0; i1 < numV; i1++) {
			InFile >> x >> y >> z;
			InFile >> col.x >> col.y >> col.z >> a;

			col /= 255;


			m_mesh->insertVertex(x,y,z);
#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
			m_mesh->vList[m_mesh->getNumV()-1].color = col;
#endif
		}

		// Read triangles
		int num_vs, v0, v1, v2, v3, v4, v5;
		for(int i1 = 0; i1 < numT; i1++) {
			InFile >> num_vs;

			//std::cerr << num_vs << std::endl;

			//if (i1%1000 == 0) std::cerr << i1 << std::endl;
			
			if (num_vs == 3) { // A triangle follows:
				InFile >> v0 >> v1 >> v2;
				//std::cerr << "Insert tri " << i1 << ": " << v0 << " " << v1 << " " << v2 << std::endl;
				m_mesh->insertTriangle(v0, v1, v2);
			} else if (num_vs == 4) {// A quadrilataral follows:
				InFile >> v0 >> v1 >> v2 >> v3;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v0, v2, v3);
			} else if (num_vs == 5) {
				InFile >> v0 >> v1 >> v2 >> v3 >> v4;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v2, v3, v4);
				m_mesh->insertTriangle(v2, v4, v0);
			} else if (num_vs == 6) {
				InFile >> v0 >> v1 >> v2 >> v3 >> v4 >> v5;
				m_mesh->insertTriangle(v0, v1, v2);
				m_mesh->insertTriangle(v2, v3, v4);
				m_mesh->insertTriangle(v2, v4, v5);
				m_mesh->insertTriangle(v2, v5, v0);
			} else if (num_vs > 6) {
				int* vs = new int[num_vs];
				for (int i1 = 0; i1 < num_vs; i1++)
					InFile >> vs[i1];

				m_mesh->insertTriangle(vs[0], vs[1], vs[2]);
				for (int i1 = 0; i1 < num_vs-4; i1++)
					m_mesh->insertTriangle(vs[2], vs[i1+3], vs[i1+4]);

				m_mesh->insertTriangle(vs[2], vs[num_vs-1], vs[0]);
				delete[] vs;
			}
		}

		//std::cerr << "Read " << numV << " Verts " << numT << " Tris" << std::endl;

		//if (numE != 0)
		//	std::cerr << "Warnung bei bool MeshReader::readOFFFile(char *filename)! Kann keine Kanten einlesen!" << std::endl;

		InFile.close();


		return m_mesh;
	};

	//******************************************************************************************

	//! Creates a SimpleMesh from an .OFF file.
	static SimpleMesh* readOFFFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();

		readOFFFile(filename, m_mesh);

		return m_mesh;
	};

	//******************************************************************************************

	//! Creates a SimpleMesh from an .OBJ file.
	static void readOBJFile(const char* filename, SimpleMesh* m_mesh) {
		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei: " << filename << std::endl;
			return;
		}

		char typ;	
		float x,y,z;
	    vec3i vertexid, textureid, normalid;
		char line[2000];
	    char tester;	

		std::vector<vec3f> normals;
		std::vector<vec2d> texturecoords;

		bool bEnd = false;
		while(!InFile.eof() || bEnd) {
			InFile >> typ;
			//fprintf(stderr, "%c  ",typ);
			if (InFile.eof()) break;
			switch(typ) {
				case 'v':
					// This may be a vertex 'v' or a vertex normal 'vn'
					if (InFile.peek() == 'n') { // dude, it's a normal
						InFile >> typ;
						InFile >> x >> y >> z;
						normals.push_back(vec3f(x,y,z));
					} else if (InFile.peek() == 't') { // dude, it's a textur coordinate
						//std::cerr << "Read textcoord" << std::endl;
						InFile >> typ;
						InFile >> x >> y;
						texturecoords.push_back(vec2d(x,y));
					} else {
						InFile >> x >> y >> z;
						m_mesh->insertVertex(x,y,z);
					}
					if (InFile.eof()) {bEnd = true; break;};
					break;
				case 'f':

					InFile >> vertexid.x;
					if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at v0\n"); break;}
					if (InFile.peek() == '/') {
						// read texture index	    
						InFile >> tester;
						if (InFile.peek() != '/') { // eskoennte sein, dass er fehlt, dass sieht dann so aus: 1 // 7
							InFile >> textureid.x; 
							if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at t0\n"); break;}
						}
					}
					if (InFile.peek() == '/') {
						// read normal index
						InFile >> tester >> normalid.x; 
						if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at n0\n"); break;}
					}

					InFile >> vertexid.y;
					if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at v1\n"); break;}
					if (InFile.peek() == '/') {
						// read texture index
						InFile >> tester;
						if (InFile.peek() != '/') {
							InFile >> textureid.y; 
							if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at t1\n"); break;}
						}
					}
					if (InFile.peek() == '/') {
						// read normal index
						InFile >> tester >> normalid.y; 
						if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at n1\n"); break;}
					}

					InFile >> vertexid.z;
					if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at v2\n"); break;}
					if (InFile.peek() == '/') {
						// read texture index
						InFile >> tester;
						if (InFile.peek() != '/') {
							InFile >> textureid.z; 
							if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at t2\n"); break;}
						}
					}
					if (InFile.peek() == '/') {
						// read normal index
						InFile >> tester >> normalid.z; 
						if (InFile.eof()) {bEnd = true; fprintf(stderr, "end at n2\n"); break;}
					}
					
					m_mesh->insertTriangle(vertexid.x-1, vertexid.y-1, vertexid.z-1);

					break;
				default:
					InFile.getline(line, 2000);	//TODO this is ugly I know
					break;
			}
		}

#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
		if (normals.size() == m_mesh->getNumV()) {
			std::cerr << "Adding normals" << std::endl;
			for (UInt i1 = 0; i1 < (UInt)m_mesh->getNumV(); i1++) {
				m_mesh->vList[i1].color = normals[i1];
			}
		}
		m_mesh->hasColors = true;
#endif

		if (texturecoords.size() == m_mesh->getNumV()) {
			std::cerr << "Adding texturecoords" << std::endl;
			for (UInt i1 = 0; i1 < (UInt)m_mesh->getNumV(); i1++) {
				m_mesh->vList[i1].param = texturecoords[i1];
			}
			m_mesh->hasTextureCoordinates = true;
		}

		InFile.close();
	}


	//! Creates a SimpleMesh from an .OBJ file.
	static SimpleMesh* readOBJFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();
	
		readOBJFile(filename, m_mesh);

		return m_mesh;
	};

//******************************************************************************************

	//! Creates a SimpleMesh from an .MESH file.
	static SimpleMesh* readMESHFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();

		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei: " << filename << std::endl;
			return NULL;
		}

		std::string typ;	
		float x,y,z;
		int t, e,f,g;
//		char line[2000];

		bool bEnd = false;
		while(!InFile.eof() || bEnd) {
			InFile >> typ;
			//fprintf(stderr, "%c  ",typ);
			if (InFile.eof()) break;

			if (typ == "Vertex") {
				InFile >> t >> x >> y >> z;
				m_mesh->insertVertex(x,y,z);
				if (InFile.eof()) {bEnd = true; break;};
			} else if (typ == "Face") {
				InFile >> t >> e >> f >> g;
				m_mesh->insertTriangle(e-1,f-1,g-1);
				if (InFile.eof()) {bEnd = true; break;};
			}
		}	

		InFile.close();
		return m_mesh;
	};

//******************************************************************************************

	//! Creates a SimpleMesh from an ascii .ASC file.
	static SimpleMesh* readGTSFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();

		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei: " << filename << std::endl;
			return NULL;
		}

		// Lese Datei Header:
		char string1[20];
		char string2[20];
		char string3[20];
		char string4[20];
		int numV, numE, numT;
		InFile >> numV >> numE >> numT >> string1 >> string2 >> string3 >> string4;


		// Read points
		double x, y, z;
		char c;
		for(int i1 = 0; i1 < numV; i1++) {
			InFile >> c >> x >> y >> z;
			m_mesh->insertVertex(x,y,z);
		}

		// Read Edges: (skip);
		int e1, e2;
		for(int i1 = 0; i1 < numE; i1++) {
			InFile >> e1 >> e2;
		}

		// Read triangles
		int v0, v1, v2;
		for(int i1 = 0; i1 < 30; i1++) {
			InFile >> v0 >> v1 >> v2;
			m_mesh->insertTriangle(v0, v1, v2);
		}

		InFile.close();
	

		return m_mesh;
	}


//******************************************************************************************

	//! Creates a SimpleMesh from an ascii .STL file.
	SimpleMesh* readSTLAsciiFile(const char* filename) {
		SimpleMesh* m_mesh = new SimpleMesh();

		std::ifstream InFile(filename);
		// Fehler abfangen
		if (InFile.fail()) {
			std::cerr << "Fehler beim Oeffnen der Datei: " << filename << std::endl;
			return NULL;
		}

		// Lese Datei Header:
		std::string currentLine;
		std::getline(InFile, currentLine);
		if (currentLine.compare("solid") != 0 && currentLine.compare("solid model") != 0) {
			std::cout << "Datei " << filename << " beginnt weder mit 'solid' noch mit 'solid model'!" << std::endl;
			return NULL;
		}

		// Lese alle Dreiecke ein:
		int line = 0;
		int triangle[3];
		for (; std::getline(InFile, currentLine) && currentLine.compare("solid") != 0; ) {
			// Entferne vorgestellte whitespaces:
			currentLine.erase(0, currentLine.find_first_not_of(" \t\n"));
			++line;

			// Aufbau einer STL ASCII Datei:
			/*
				"facet normal X X X"
				"outer loop"
				"vertex X X X"
				"vertex X X X"
				"vertex X X X"
				"endloop"
				"endfacet"
			*/
			if (   (currentLine.compare(0, 12, "facet normal", 0, 12) == 0)
				|| (currentLine.compare(0, 10, "outer loop", 0, 10) == 0)
				|| (currentLine.compare(0, 7, "endloop", 0, 7) == 0)) {
				// Ignoriere Zeile
				continue;
			}

			// Dreieck einlesen abgeschlossen; sechste Zeile lautet:
			if (currentLine.compare(0, 8, "endfacet", 0, 8) == 0) {
				m_mesh->insertTriangle(triangle[0], triangle[1], triangle[2]);
				line = 0;
				continue;
			}

			// Zeile 3-5 lauten:
			// "vertex X X X"
			// Zerlege String in doubles
			double x, y, z;
			std::string tmpString;
			std::istringstream istream;
			istream.clear();
			istream.str(currentLine);
			istream >> tmpString >> x >> y >> z;
			// Füge Punkt mit den entsprechenden Koordinaten in das Mesh ein:
			m_mesh->insertVertex(x, y, z);
			triangle[line-3] = m_mesh->getNumV()-1;
		}
		InFile.close();

		return m_mesh;
	}

protected:

private:
};


