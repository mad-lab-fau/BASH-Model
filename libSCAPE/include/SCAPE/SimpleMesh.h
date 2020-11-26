#pragma once

typedef unsigned int UInt;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 //declare our M_PI constant
#endif

#ifndef M_PI2
#define M_PI2 6.2831853 //declare our M_PI2 constant
#endif

#include <vector>
#include <set>

#include "point3d.h"
#include "point2d.h"


class Triangle;
class TriangleList;
class Edge;
class EdgeList;
class Vertex;
class SimpleMesh;


void UmbrellaSmooth(SimpleMesh* sm);

//******************************************************************************************
//************************************** OBB ***********************************************
//******************************************************************************************

class edgeCollapsInfo {

public:

	//! sm->vlist[basis[0]].c, sm->vlist[basis[1]].c and sm->vlist[basis[2]].c span the basis for the newly added point
	vec3i basis;
	//! local coordinates to the 'to be added' point in the above basis.
	vec3d local_coords;
	//! The vertices of the fan in which the new vertex will be inserted.
	std::vector<vec2i> fan_edges;
	//! tri_corner_0 + two succeeding vertices of 'triangleIDstoRemove' mark the triangles which have been filled into the hole.
	int tri_corner_0;

	//! marks whether the removed vertex was a boundary vertex
	bool isBoundaryVertex;

};

//******************************************************************************************
//************************************** OBB ***********************************************
//******************************************************************************************


//! an object aligned Bounding Box.
template<class T>
class OBB {

public:

	OBB() {}

	~OBB() {}

	//! Center of the OBB.
	point3d<T> center;
	//! Basis of the OBB.
	T basis[9];
	//! Basis Vector in x-dir of the OBB.
	point3d<T> base_x;
	//! Basis Vector in y-dir of the OBB.
	point3d<T> base_y;
	//! Basis Vector in z-dir of the OBB.
	point3d<T> base_z;
	//! Extend in x-dir of basis
	T dx;
	//! Extend in y-dir of basis
	T dy;
	//! Extend in z-dir of basis
	T dz;

	void print() {
		std::cerr << "dx: " << dx << std::endl;
		std::cerr << "dy: " << dy << std::endl;
		std::cerr << "dz: " << dz << std::endl;
	}

	//! Copies values of base_x, base_y and base_x to basis
	void copyToBasis() {
		basis[0] = base_x.x;
		basis[1] = base_y.x;
		basis[2] = base_z.x;
		basis[3] = base_x.y;
		basis[4] = base_y.y;
		basis[5] = base_z.y;
		basis[6] = base_x.z;
		basis[7] = base_y.z;
		basis[8] = base_z.z;
	}

	//! Copies values of base_x, base_y and base_x to scaled basis
	void copyToScaledBasis() {
		basis[0] = base_x.x*dx;
		basis[1] = base_y.x*dy;
		basis[2] = base_z.x*dz;
		basis[3] = base_x.y*dx;
		basis[4] = base_y.y*dy;
		basis[5] = base_z.y*dz;
		basis[6] = base_x.z*dx;
		basis[7] = base_y.z*dy;
		basis[8] = base_z.z*dz;
	}

};


//******************************************************************************************
//*********************************** TRIANGLE *********************************************
//******************************************************************************************


//! a triangle that will be used for the Mesh
class Triangle
{
public:

	friend class SimpleMesh;

	// Constructor.
	Triangle(void);

	// Constructor.
	Triangle(int v0, int v1, int v2, Edge *e0, Edge *e1, Edge *e2, SimpleMesh* mesh);

	// Destructor.
	~Triangle(void);

	//! returns m_v0. 
	inline int v0() const { return m_v0;	};
	//! returns m_v1. 
	inline int v1() const { return m_v1;	};
	//! returns m_v2. 
	inline int v2() const { return m_v2;	};

	inline int getV( int x ) const {
		if (x == 0)
			return m_v0;
		if (x == 1)
			return m_v1;
		// else
		return m_v2;
	}

	//! Returns the interpolated color using the barycentric coordinates barys.
	vec3f getInterpolatedColor(const vec3d& barys) const;

	//! Returns the interpolated pos using the barycentric coordinates barys.
	vec3d getInterpolatedPos(const vec3d& barys) const;

	//! Computes and returns the Normal of the Triangle.
	vec3d getNormal() const;

	//! Computes and returns the Area of the Triangle.
	double getArea() const;

	//! returns the barycentric coordinates of the point "p" w.r.t this triangle (assuming the triangle is 2D, i.e. z=0 for all vertices).
	vec3d barycentric(vec2d p) const;

	//! computes and returns the barycentric coordinates of the point "p" w.r.t this triangle after projecting the point onto the tangentplane of this triangle
	vec3d barycentric(vec3d p) const;

	//! [WITH SIDE EFFECTS] computes and returns the barycentric coordinates of the point "p" w.r.t this triangle after projecting the point onto the tangentplane of this triangle
	/*
		This function <b>changes 'p'<b> by projecting 'p' onto the tangentplane of the current triangle!!
	*/
	vec3d barycentric_with_proj(vec3d& p) const;

	//! returns the center-point of this triangle 1/3(v0+v1+v2);
	vec3d getCenter() const;

	//! Access v0, v1, v2.
	int& operator[](unsigned int i);

	//! Returns the number of corresponding Edges.
	/*!
		compare == 1 means the two triangles share a vertex.<br>
		compare == 2 means the two triangles share an edge.<br>
		compare == 3 means the two triangles are equal.<br>
	*/
	int compare(Triangle other);

	//! Returns the i'th (0...2) neighbortriangle or NULL;
	Triangle* getTriangleNeighbor(int i);

	//! Returns the third vertex of the triangle.
	int getOther(int v0, int v1);

	//! If the vertex 'v' is contained, the angle at 'v' in the triangle is returned, else 0.
	double getAngleAtVertex(int v);

	//! Prints the vertices of the triangle to stderr.
	void print();

	//! TODO: Why are these Edges public?
	Edge *e0;
	//! TODO: Why are these Edges public?
	Edge *e1;
	//! TODO: Why are these Edges public?
	Edge *e2;

	//! A marker that can be used for various tasks (For example in 'operators/meshMadgicWand')
	bool marker;

#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
	int intFlag;
#endif

	//! Changes the orientation of the Triangle from CW to CCW and vice-versa.
	void invertTriangle() {
		int v0 = m_v0;
		int v1 = m_v1;
		m_v0 = v1;
		m_v1 = v0;
	};

	void set_v0(int v) {
		m_v0 = v;
	};
	void set_v1(int v) {
		m_v1 = v;
	};
	void set_v2(int v) {
		m_v2 = v;
	};

	//! Returns if the vertices of the edges match the vertices of the triangle.
	bool isConsistent() const;

	//! The indicees of the points of the triangles as stored in m_mesh.
	int m_v0, m_v1, m_v2;
private:
	//! The parent-mesh
	SimpleMesh* m_mesh;

};


//******************************************************************************************
//******************************** TRIANGLELIST ********************************************
//******************************************************************************************


//! vector of Triangle*.
class TriangleList : public std::vector<Triangle*>
{
public:

	//! Destructor
	~TriangleList() {
		//for (UInt i1 = 0; i1 < size(); i1++) delete (*this)[i1];
	};

	//! inserts a triangle, i.e. calls push_back(t).
	void insertTriangle(Triangle* t);

	//! removes the triangle(v0, v1, v2) from the list.
	/*!
		\return Returns whether the triangle has been removed.
	*/
	bool removeTriangle(int v0, int v1, int v2);

	//! returns a pointer to the triangle(v0, v1, v2) if it's the list.
	/*!
		\return The triangle or NULL.
	*/
	Triangle* getTriangle(int v0, int v1, int v2);

	//! returns the index of the triangle(v0, v1, v2) if it's the list.
	/*!
		\return Index of the triangle of 0
	*/
	UInt getTriangleIndex(int v0, int v1, int v2);

	//! Prints to stdcerr.
	void print() {
		for (UInt i1 = 0; i1 < size(); i1++) operator[](i1)->print();
	}
};


//******************************************************************************************
//************************************* EDGE ***********************************************
//******************************************************************************************


//! edge of a mesh
class Edge
{
public:

	//! Constructor.
	Edge(void);

	//! Constructor.
	Edge(int v0_, int v1_, SimpleMesh* mesh);

	//! Destructor.
	~Edge(void);

	//! Returns true if the two Edges have the same endpoints.
	bool operator==(const Edge& other);

	//! Returns the other vertex of the edge.
	inline int getOther(int v) const {
		return ((v == v0) ? v1 : v0);
	}

	//! Returns whether the edge is convex (true) or concave (false).
	bool isConvex() {
		if (tList.size() != 2) {
			std::cerr << "It doesn't make sense to call Edge::isConvex() on an Edge with " << (UInt)tList.size() << " adjacent triangles" << std::endl;
			return false;
		}

		Triangle* t0 = tList[0];
		Triangle* t1 = tList[1];
		vec3d n0 = t0->getNormal();
		vec3d n1 = t1->getNormal();
		vec3d c0 = t0->getCenter();
		vec3d c1 = t1->getCenter();

		vec3d c0c1 = c1-c0;
		double angle0 = c0c1|n0;
		double angle1 = -c0c1|n1;
		double a0 = acos(angle0);
		double a1 = acos(angle1);

		if (a0 + a1 <= M_PI) return true;
		else return false;

	}

	//! Returns the index of the 'i'th triangle in the array m_mesh->tList.
	UInt getIndexOfTriangle(UInt i);

	//! Swaps the edge (If this decreases the curvature of the mesh) and tells neighbors to do the same recursively.
	void Swap();

	//! Prints endpoints to stderr.
	void print();

	//! returns the length of the edge.
	double getLength() const;

	//! Computes the sum of the dihedral angles at the edges and the edges of the adjacent triangles.
	double computeDihedralAngles() const;

	//! Computes the sum of the dihedral angle at the edge.
	double computeDihedralAngle() const;

	//! Computes the cotangent weights at this edge.
	void computeCotangentWeights(bool uniform = false);

	void computeCotangentWeightsVERBOSE() const;

	//! Can be used to indicate for example that the egde is sharp.
	bool marker;

	//! The triangles adjacing to this edge.
	TriangleList tList;

	//! The parent-mesh
	SimpleMesh* m_mesh;

	//! The indicees of the points of the triangles as stored in m_mesh.
	int v0, v1;

	//! the cotangent-weights
	double cotangent_weight;
};

//******************************************************************************************
//*********************************** EDGELIST *********************************************
//******************************************************************************************


//! vector of edges
class EdgeList : public std::vector<Edge*>
{
public:

	~EdgeList() {
		//for (UInt i1 = 0; i1 < size(); i1++) delete (*this)[i1];
	};

	//! Returns the edge with endpoints v0 and v1 if it is in the list, else NULL.
	Edge* getEdge(int v0, int v1);

	//! Inserts an edge, i.e. calls push_back(e).
	void insertEdge(Edge* e);

	//! Prints the list.
	void print();

	//! removes the edge(v0, v1) from the list.
	/*!
		\return Returns whether the edge has been removed.
	*/
	bool removeEdge(int v0, int v1);

	//! removes the edge 'e' from the list.
	/*!
		\return Returns whether the edge has been removed.
	*/
	bool removeEdge(Edge* e);

	//! The parent-mesh.
	SimpleMesh* m_mesh;

};


//******************************************************************************************
//************************************* VERTEX *********************************************
//******************************************************************************************


//! Vertex of a mesh.
class Vertex
{
public:

	//! Constructor
	Vertex(void);

	//! Constructor
	Vertex(vec3d c_);

	//! Destructor
	~Vertex(void);

	//! The coordinates (x,y,z) of the vertex.
	vec3d c;

	//! The parametrization of the vertex.
	vec2d param;

	//! All edges with this vertex as endpoint.
	EdgeList eList;

#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
	//! The normal of the vertex <b>Caution, same as normal</b>.
	vec3f color;
	//union {
	//	//! The color of the vertex <b>Caution, same as color</b>.
	//	vec3f normal;
	//};
#endif

	//! Returns the number of neighbors of the vertex (which is eList.size()).
	inline UInt getNumN() const {
		return (UInt) eList.size();
	};

	//! Returns the 'j'-th neighbor of the vertex '#'.
	inline int getNeighbor(UInt j, int i) const {
		return eList[j]->getOther(i);
	}

	//! Returns a set containing all triangles adjacent to the vertex.
	std::set<Triangle*> getAdjacentTriangles() {

		std::set<Triangle*> tList;
		for (UInt i1 = 0; i1 < getNumN(); i1++) {
			for (UInt i2 = 0; i2 < eList[i1]->tList.size(); i2++) {
				Triangle* t = eList[i1]->tList[i2];
				tList.insert(t);
			}
		}

		return tList;
	}

	union {
		//! Marks whether the vertex lies on the boundary of the mesh. <b>Caution, same as has_trace and bool_flag</b> Initialized to 'false'.
		bool is_boundary_point;
		//! Marks whether the vertex still has trace (is used for deformable registration). <b>Caution, same as is_boundary_point and bool_flag</b> Initialized to 'false'.
		bool has_trace;
		//! A boolean flag which may be used to assign a boolean value to a vertex. <b>Caution, same as is_boundary_point and has_trace</b>. Initialized to 'false'.
		bool bool_flag;
	};
};


//******************************************************************************************
//********************************** SIMPLE MESH *******************************************
//******************************************************************************************


//! A (not so) simple mesh class
class SimpleMesh
{

	//! Allow the function closeSimpleHoles(...) to access our private member functions
	friend void closeSimpleHoles(SimpleMesh* sm);


public:
	//! Constructor.
	SimpleMesh(void);

	//! Copy Constructor
	SimpleMesh(const SimpleMesh& other) {
		for (int i1 = 0; i1 < other.getNumV(); i1++) {
			insertVertex(other.vList[i1].c);
		}
		for (int i1 = 0; i1 < other.getNumT(); i1++) {
			insertTriangle(other.tList[i1]->m_v0, other.tList[i1]->m_v1, other.tList[i1]->m_v2);
		}
		setNumV(other.getNumV());
		m_numT = other.m_numT;
		hasTextureCoordinates = false;
	}

	//! Copy Assignment Operator
	SimpleMesh& operator=(const SimpleMesh& other) {
		vList.clear();
		tList.clear();
		eList.clear();
		for (int i1 = 0; i1 < other.getNumV(); i1++) {
			insertVertex(other.vList[i1].c);
		}
		for (int i1 = 0; i1 < other.getNumT(); i1++) {
			insertTriangle(other.tList[i1]->m_v0, other.tList[i1]->m_v1, other.tList[i1]->m_v2);
		}
		setNumV(other.getNumV());
		m_numT = other.m_numT;
		hasTextureCoordinates = false;

		return *this;
	}

	SimpleMesh(const float* vertices, const UInt* indices, int numVertices, int numIndices) {
		m_numV = 0;
		m_numT = 0;
		boundary_vertices_computed = false;
		hasTextureCoordinates = false;

		for (int i = 0; i < numVertices * 3; i += 3) {
			float x = vertices[i];
			float y = vertices[i+1];
			float z = vertices[i+2];
			insertVertex(x, y, z);
		}
		for (int i = 0; i < numIndices; i += 3) {
			UInt i0 = indices[i];
			UInt i1 = indices[i+1];
			UInt i2 = indices[i+2];
			insertTriangle(i0, i1, i2);
		}
	}

	template<class V, class T>
	SimpleMesh(const std::vector< V >& pos, const std::vector< T >& tris) {
		m_numV = 0;
		m_numT = 0;
		boundary_vertices_computed = false;
		hasTextureCoordinates = false;

		for (UInt i1 = 0; i1 < pos.size(); i1++) insertVertex(pos[i1].x, pos[i1].y, pos[i1].z);
		for (UInt i1 = 0; i1 < tris.size(); i1++) insertTriangle(tris[i1].x, tris[i1].y, tris[i1].z);
	}

	//! Destructor.
	~SimpleMesh(void);

	//! Reads an .OFF file.
	//bool readOffFile(char* filename);

	//! Inserts a Edge with the given indicees.
	Edge* insertEdge(int v0, int v1);

	//! Inserts a Triangle with the given indicees, does also update neighbor information.
	Triangle* insertTriangle(int v0, int v1, int v2);

	//! returns a pointer to the triangle(v0, v1, v2) if it's contained in the mesh.
	/*!
		\return The triangle or NULL.
	*/
	inline Triangle* getTriangle(int v0, int v1, int v2) {
		Edge* e = vList[v0].eList.getEdge(v0, v1);
		if (e == NULL) return NULL;
		return (e->tList.getTriangle(v0,v1,v2));
	};

	//! Inserts a Vertex.
	void insertVertex(double x, double y, double z);

	//! Inserts a Vertex.
	inline void insertVertex(vec3d vert) {
		insertVertex(vert.x, vert.y, vert.z);
	};

	//! Returns the coordinates of the "i"th vertex.
	inline vec3d getVertex(int i) const {
		if (i >= m_numV) {
			std::cerr << "Out of bounds access on m_points" << std::endl;
			exit(EXIT_FAILURE);
		}
		return vList[i].c;
	}

	//! Prints the whole mesh to stderr.
	void print();

	//! Returns the boundary of the mesh.
	std::vector<int> getBoundary();

	//! Returns the boundary of the mesh starting at edge 'e'.
	std::vector<int> getBoundaryStartAtEdge(Edge* e, bool& is_included);

	//! Returns the number of Triangles in the mesh.
	inline int getNumT() const {
		return m_numT;
	};

	//! Returns the number of Vertices in the mesh.
	inline int getNumV() const {
		return (int)vList.size();
	};

	//! Returns the number of edges that belong to only one triangle
	inline int getNumOpenEdges() const {
		int numOpenEdges = 0;
		for (unsigned int i1 = 0; i1 < eList.size(); i1++) {
			if (eList[i1]->tList.size() < 2)
                numOpenEdges++;
		}
		return numOpenEdges;
	};

	//! Returns the number of edges which do not have 2 neighbors (i.e. 0, 1, 3, 4, 5, ..., neighbors)
	inline int getNumBadEdges() const {
		int numBadEdges = 0;
		for (unsigned int i1 = 0; i1 < eList.size(); i1++) {
			if (eList[i1]->tList.size() != 2) {
                numBadEdges++;
				//if (eList[i1]->tList.size() > 2)
				//	std::cerr << (UInt)eList[i1]->tList.size() << " neighbors at edge " << eList[i1]->v0 << " / " << eList[i1]->v1 << std::endl;
			}
		}
		return numBadEdges;
	};

	//! Recomputes all vertices as pt_new = (pt+move)*scale;
	inline void resize(double scale, vec3d move) {
		for (int i1 = 0; i1 < getNumV(); i1++) {
			vList[i1].c = (vList[i1].c+move)*scale;
		}
	}

	//! Manipulate private member m_numV
	void setNumV(int n) {
		m_numV = n;
	}

	//! Marks all boundaryvertices.
	void computeBoundaryVertices();

	//! Removes the triangle 't' from the list. Does also take care about the edges and vertices.
	void removeTriangle(Triangle* t);

	//! removes the triangle(v0, v1, v2) from the eLists of its edges. Does not delete the triangle physically but sets t->marker=true
	void removeTriangleNoPhysical(Triangle* t);

	//! removes the zombies of 'removeTriangleNoPhysical' from the global tList.
	void cleanUpMarkedTriangles();

	//! Removes the Vertex 'id' and all it's edges and triangles.
	void removeVertex(int id);

	//! Removes the Vertex 'id' and all it's edges and triangles. Calls removeVertexPhysically() to remove the vertex from the vList.
	void removeVertexAndPhysically(int id) {
		removeVertex(id);
		removeVertexPhysically(id);
	}

	//! Removes all disconnected vertices.
	void cleanValence0Vertices(std::vector<edgeCollapsInfo>* eColInfo = NULL);

	//! Removes the edge 'v0'-'v1' and all adjacent triangles.
	void removeEdge(int v0, int v1);

	//! Removes the edge 'e' and all adjacent triangles.
	void removeEdge(Edge* e);

	//! Removes the vertex 'v0' by performing the optimal edge-collaps.
	void edgeCollapse(int v0, edgeCollapsInfo* colInfo = NULL);

	//! Computes the Object-Aligned Bounding Box of the mesh.
	void computeOBB();

	//! returns the centroid of the fan around vertex 'v'.
	vec3d getFanCentroid(int v);

	//! Recomputes 'min' and 'max'.
	void recomputeMinMax() {
		min = vList[0].c;
		max = vList[0].c;

		for (UInt i1 = 1; i1 < vList.size(); i1++) {
			const vec3d& pt = vList[i1].c;
			if(pt.x < min.x) min.x = pt.x;
			if(pt.x > max.x) max.x = pt.x;
			if(pt.y < min.y) min.y = pt.y;
			if(pt.y > max.y) max.y = pt.y;
			if(pt.z < min.z) min.z = pt.z;
			if(pt.z > max.z) max.z = pt.z;
		}
	}


//******************************************************************************************
//	Member Variables
//******************************************************************************************


	//! The vertices of the mesh.
	std::vector<Vertex> vList;

	//! The triangles of the mesh.
	TriangleList tList;

	//! The edges of the mesh.
	EdgeList eList;

	//! Extend of the mesh, min[0] denotes the smallest x value, max[1] the largest y value...
	vec3d min, max;

	//! The BoundingBox of the Mesh, need to be explicitely computed by calling 'computeOBB'.
	OBB<double> boundingbox;

	//! set if the mesh has texture coordinates
	bool hasTextureCoordinates;

	//! may be set if the mesh has colors
	bool hasColors;
protected:

	//! Removes 'id' from 'vList', decreases all references to vertices with an id > 'id'.
	void removeVertexPhysically(int id);

	//! Number of Triangles
	int m_numT;
	//! Number of Points
	int m_numV;

	//! Will be set to true after the boundary vertices have been computed.
	bool boundary_vertices_computed;
};

