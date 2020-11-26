#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "../SimpleMesh.h"


//****************************************************************************************************


//! Writes a SimpleMesh to a file.
class MeshWriter
{
public:


//****************************************************************************************************

	//! Constructor.
	MeshWriter(){};

//****************************************************************************************************

	//! Writes an .OFF file from a a vector of indices and triangles.
	template <class VERTEXTYPE, class FLOATTYPE>
	static void writeOFFFile(const char* filename, const std::vector<VERTEXTYPE>& verts, const std::vector<std::vector<int> >& tris, FLOATTYPE scale = 1, point3d<FLOATTYPE> move = point3d<FLOATTYPE>(0,0,0)) {
		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "OFF" << std::endl;
		int numT = (int)tris.size();
		int numV = (int) verts.size();
		int numE = 0;

		OutFile << numV << " " << numT <<  " " << numE << std::endl;

		for (size_t i1 = 0; i1 < verts.size(); i1++) {
			const VERTEXTYPE& pt = verts[i1];
			OutFile << (pt.point[0]+move[0])*scale << " " << (pt.point[1]+move[1])*scale << " " << (pt.point[2]+move[2])*scale << std::endl;
		}

		for (size_t i1 = 0; i1 < tris.size(); i1++) {
			const std::vector<int>& tri = tris[i1];
			OutFile << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << std::endl;
		}
}

//****************************************************************************************************

	//! Writes an .OBJ file from a SimpleMesh, uses the information stored in the color field as normals.
	static void writeOBJFileWithNormals(const char* filename, SimpleMesh* m_mesh) {
		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "g points" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();

		//std::cerr << "Write : " << numV << " Vertices and " << numT << " Tris" << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << "v " << m_mesh->vList[i1].c.x <<  " " << m_mesh->vList[i1].c.y <<  " " <<  m_mesh->vList[i1].c.z <<  std::endl;
		}

#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
		// write normals
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << "vn " << m_mesh->vList[i1].color.x <<  " " << m_mesh->vList[i1].color.y <<  " " <<  m_mesh->vList[i1].color.z <<  std::endl;
		}
#endif

		OutFile << "g surface" << std::endl;

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
			OutFile <<  "f " << m_mesh->tList[i1]->v0()+1 << " " << m_mesh->tList[i1]->v1()+1 << " " << m_mesh->tList[i1]->v2()+1 << std::endl;
		}

		OutFile.close();
		//std::cerr << "Wrote file " << filename << std::endl;
	}

//****************************************************************************************************

	//! Writes an .OBJ file from a SimpleMesh.
	static void writeOBJFileWithTextureCoords(const char* filename, SimpleMesh* m_mesh) {

		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "g points" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();

		//std::cerr << "Write : " << numV << " Vertices and " << numT << " Tris" << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << "v " << m_mesh->vList[i1].c.x <<  " " << m_mesh->vList[i1].c.y <<  " " <<  m_mesh->vList[i1].c.z <<  std::endl;
		}
		// write texture coordinates
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << "vt " << m_mesh->vList[i1].param.x <<  " " << m_mesh->vList[i1].param.y <<  std::endl;
		}
		OutFile << "g surface" << std::endl;

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
			//OutFile <<  "f " << m_mesh->tList[i1]->v0()+1 << " " << m_mesh->tList[i1]->v1()+1 << " " << m_mesh->tList[i1]->v2()+1 << std::endl;
			OutFile <<  "f " << m_mesh->tList[i1]->v0()+1 << "/"  << m_mesh->tList[i1]->v0()+1 << " " << m_mesh->tList[i1]->v1()+1 << "/"  << m_mesh->tList[i1]->v1()+1  << " " << m_mesh->tList[i1]->v2()+1 << "/"  << m_mesh->tList[i1]->v2()+1  << std::endl;
		}

		OutFile.close();
		//std::cerr << "Wrote file " << filename << std::endl;

	};

//****************************************************************************************************

	//! Writes an .OBJ file from a SimpleMesh.
	static void writeOBJFile(const char* filename, SimpleMesh* m_mesh) {

		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "g points" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();

		//std::cerr << "Write : " << numV << " Vertices and " << numT << " Tris" << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << "v " << m_mesh->vList[i1].c.x <<  " " << m_mesh->vList[i1].c.y <<  " " <<  m_mesh->vList[i1].c.z <<  std::endl;
		}

		OutFile << "g surface" << std::endl;

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
			OutFile <<  "f " << m_mesh->tList[i1]->v0()+1 << " " << m_mesh->tList[i1]->v1()+1 << " " << m_mesh->tList[i1]->v2()+1 << std::endl;
		}

		OutFile.close();
		//std::cerr << "Wrote file " << filename << std::endl;

	};

	//****************************************************************************************************

	//! Writes an .BIN file from a simple mesh.
	static void writeBINFile(const char* filename, SimpleMesh* m_mesh, bool writeOnlyMarkedTriangles=false, bool markedWithTrue=true) {

		std::ofstream OutFile(filename, std::ios::binary);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();
		int numT_toWrite = numT;
		if (writeOnlyMarkedTriangles) {
			numT_toWrite = 0;
			for(int i1 = 0; i1 < numT; i1++) {
				if (m_mesh->tList[i1]->marker == markedWithTrue)
					numT_toWrite++;
			}
		}
		int numE = 0;

		OutFile.write((char*)&numV, sizeof(int));
		OutFile.write((char*)&numT_toWrite, sizeof(int));

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			double x = m_mesh->vList[i1].c[0];
			double y = m_mesh->vList[i1].c[1];
			double z = m_mesh->vList[i1].c[2];
			// catch NANs ind INFs.
			if (x != x || x == std::numeric_limits<double>::infinity() || x == -std::numeric_limits<double>::infinity() ) x = 0;
			if (y != y || y == std::numeric_limits<double>::infinity() || y == -std::numeric_limits<double>::infinity() ) y = 0;
			if (z != z || z == std::numeric_limits<double>::infinity() || z == -std::numeric_limits<double>::infinity() ) z = 0;

			OutFile.write((char*)&m_mesh->vList[i1].c[0], 3*sizeof(double));
		}

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
			if (!writeOnlyMarkedTriangles || (m_mesh->tList[i1]->marker==markedWithTrue)) {
				vec3i t(m_mesh->tList[i1]->v0(), m_mesh->tList[i1]->v1(), m_mesh->tList[i1]->v2());
				OutFile.write((char*)&t, 3*sizeof(int));
			}
		}


		OutFile.close();
		std::cerr << "Wrote file " << filename << std::endl;

	};

//****************************************************************************************************

	//! Writes an .OFF file from a simple mesh with color.
	static void writeOFFFileVertexColor(const char* filename, SimpleMesh* m_mesh) {

		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "COFF" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();
		int numT_toWrite = numT;
		int numE = 0;

		OutFile << numV << " " << numT_toWrite <<  " " << numE << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			double x = m_mesh->vList[i1].c[0];
			double y = m_mesh->vList[i1].c[1];
			double z = m_mesh->vList[i1].c[2];

			OutFile << x <<  " " << y <<  " " << z;
			vec3f c = m_mesh->vList[i1].color;
			c *= 255;
			OutFile << " " << c.x <<  " " << c.y <<  " " << c.z << " 255" << std::endl;
		}

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			OutFile <<  "3 " << m_mesh->tList[i1]->v0() << " " << m_mesh->tList[i1]->v1() << " " << m_mesh->tList[i1]->v2() << std::endl;
		}


		OutFile.close();
		std::cerr << "Wrote file " << filename << std::endl;
	}

//****************************************************************************************************

	//! Writes an .OFF file from a simple mesh.
	static void writeOFFFile(const char* filename, SimpleMesh* m_mesh, bool writeOnlyMarkedTriangles=false, bool markedWithTrue=true) {

		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "OFF" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();
		int numT_toWrite = numT;
		if (writeOnlyMarkedTriangles) {
			numT_toWrite = 0;
			for(int i1 = 0; i1 < numT; i1++) {
				if (m_mesh->tList[i1]->marker == markedWithTrue)
					numT_toWrite++;
			}
		}
		int numE = 0;

		OutFile << numV << " " << numT_toWrite <<  " " << numE << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			double x = m_mesh->vList[i1].c[0];
			double y = m_mesh->vList[i1].c[1];
			double z = m_mesh->vList[i1].c[2];
			// catch NANs ind INFs.
			if (x != x || x == std::numeric_limits<double>::infinity() || x == -std::numeric_limits<double>::infinity() ) x = 0;
			if (y != y || y == std::numeric_limits<double>::infinity() || y == -std::numeric_limits<double>::infinity() ) y = 0;
			if (z != z || z == std::numeric_limits<double>::infinity() || z == -std::numeric_limits<double>::infinity() ) z = 0;

			OutFile << x <<  " " << y <<  " " << z <<  std::endl;
		}

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
			if (!writeOnlyMarkedTriangles || (m_mesh->tList[i1]->marker==markedWithTrue))
				OutFile <<  "3 " << m_mesh->tList[i1]->v0() << " " << m_mesh->tList[i1]->v1() << " " << m_mesh->tList[i1]->v2() << std::endl;
		}


		OutFile.close();
		std::cerr << "Wrote file " << filename << std::endl;

	};

//****************************************************************************************************

	//! Writes an .OFF file from a simple mesh.
	static void writeOFFFileChart(const char* filename, SimpleMesh* m_mesh, int chartID) {

		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			std::cerr << "Could not write to [" << filename << "]" << std::endl;
			return;
		}

		OutFile << "OFF" << std::endl;

		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();
		int numT_toWrite = 0;
#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
		for(int i1 = 0; i1 < numT; i1++) {
			if (m_mesh->tList[i1]->intFlag == chartID)
				numT_toWrite++;
		}
#endif
		int numE = 0;

		OutFile << numV << " " << numT_toWrite <<  " " << numE << std::endl;

		// write points
		for(int i1 = 0; i1 < numV; i1++) {
			double x = m_mesh->vList[i1].c[0];
			double y = m_mesh->vList[i1].c[1];
			double z = m_mesh->vList[i1].c[2];
			// catch NANs ind INFs.
			if (x != x || x == std::numeric_limits<double>::infinity() || x == -std::numeric_limits<double>::infinity() ) x = 0;
			if (y != y || y == std::numeric_limits<double>::infinity() || y == -std::numeric_limits<double>::infinity() ) y = 0;
			if (z != z || z == std::numeric_limits<double>::infinity() || z == -std::numeric_limits<double>::infinity() ) z = 0;

			OutFile << x <<  " " << y <<  " " << z <<  std::endl;
		}

		// write triangles
		for(int i1 = 0; i1 < numT; i1++) {
			// std::cerr << "Writing Triangle " << i1 << " of " << numT << std::endl;
#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
	if (m_mesh->tList[i1]->intFlag == chartID)
#endif
				OutFile <<  "3 " << m_mesh->tList[i1]->v0() << " " << m_mesh->tList[i1]->v1() << " " << m_mesh->tList[i1]->v2() << std::endl;
		}


		OutFile.close();
		std::cerr << "Wrote file " << filename << std::endl;

	};

//****************************************************************************************************

	//! Writes an binary .STL file from a simple mesh.
	static void writeSTLBinaryFile(const char* filename, SimpleMesh* sm) {
		std::ofstream OutFile(filename, std::ios::binary);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			return;
		}

		// Write Header
		OutFile << "STL binary by JBS_lib                   ";
		OutFile << "                              end Header";
        
		UInt numT = sm->getNumT();
		OutFile.write((char*)&numT, sizeof(unsigned int));

		short x = 0;

		for (UInt i1 = 0; i1 < numT; i1++) {
			vec3d normal = sm->tList[i1]->getNormal();
			vec3d v0 = sm->vList[sm->tList[i1]->v0()].c;
			vec3d v1 = sm->vList[sm->tList[i1]->v1()].c;
			vec3d v2 = sm->vList[sm->tList[i1]->v2()].c;

			vec3f f_normal((float)normal.x, (float)normal.y, (float)normal.z);
			vec3f f_v0((float)v0.x, (float)v0.y, (float)v0.z);
			vec3f f_v1((float)v1.x, (float)v1.y, (float)v1.z);
			vec3f f_v2((float)v2.x, (float)v2.y, (float)v2.z);

			OutFile.write((char*)&f_normal, 3*sizeof(float));
			OutFile.write((char*)&f_v0, 3*sizeof(float));
			OutFile.write((char*)&f_v1, 3*sizeof(float));
			OutFile.write((char*)&f_v2, 3*sizeof(float));

			OutFile.write((char*)&x, sizeof(short));
		}

		OutFile.close();

	}

//****************************************************************************************************

//****************************************************************************************************

	//! Writes an .IV file from a simple mesh.
	static void writeIVFile(const char* filename, SimpleMesh* m_mesh) {
		std::ofstream OutFile(filename);
		// Fehler abfangen
		if (OutFile.fail()) {
			std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
			return;
		}

		OutFile << "#Inventor V2.1 ascii" << std::endl;
		OutFile << std::endl;
		OutFile << "Separator {" << std::endl;
		OutFile << "	Coordinate3 {" << std::endl;
		OutFile << "		point [" << std::endl;
		int numV = m_mesh->getNumV();
		int numT = m_mesh->getNumT();
		for(int i1 = 0; i1 < numV; i1++) {
			OutFile << m_mesh->vList[i1].c[0] <<  " " << m_mesh->vList[i1].c[1] <<  " " <<  m_mesh->vList[i1].c[2] << "," << std::endl;
		}
		OutFile << "		]" << std::endl;
		OutFile << "	}" << std::endl;
		OutFile << "IndexedFaceSet {" << std::endl;
		OutFile << "		coordIndex[" << std::endl;
		for(int i1 = 0; i1 < numT; i1++) {
			OutFile << m_mesh->tList[i1]->v0() << ", " << m_mesh->tList[i1]->v1() << ", " << m_mesh->tList[i1]->v2() << ", -1," << std::endl;
		}
		OutFile << "		]" << std::endl;
		OutFile << "	}" << std::endl;
		OutFile << "}" << std::endl;

		OutFile.close();
	}

//****************************************************************************************************


protected:

private:
};


