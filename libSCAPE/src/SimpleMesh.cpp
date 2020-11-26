/*
   An explanation on the "MISSING:" and "DUE TO MISSING:" comments
   ================================================================

   Some include files are missing, therefore referenced functions will not be found.
   As a first step to get the rest of this code working, functions in this this code unit
   which references unknown identifiers, possibly from these non-existing include files,
   are generously commented out.
*/

#include "SimpleMesh.h"
#include "MatrixnD.h"

#include <set>


//******************************************************************************************
//*********************************** ******** *********************************************
//******************************************************************************************

struct TripleInt {
	TripleInt(int v0_, int v1_, int v2_) : v0(v0_), v1(v1_), v2(v2_) {}
	int v0, v1, v2;
};

//******************************************************************************************
//*********************************** ******** *********************************************
//******************************************************************************************


int getMinimumTriangulation(std::vector< std::vector<TripleInt> >& triangulations, SimpleMesh* sm) {
	int		best_id = 0;
	double	best_weight = std::numeric_limits<double>::max();
	//std::cerr << "triangulations.size(): " << triangulations.size() << std::endl;
	//std::cerr << "numTri: " << triangulations[0].size() << std::endl;
	for (UInt i1 = 0; i1 < triangulations.size(); i1++) {
		const std::vector<TripleInt>& tri = triangulations[i1];
		double w = 0;
		for (UInt i2 = 0; i2 < tri.size()-1; i2++) {
			const TripleInt& t1 = tri[i2];
			const TripleInt& t2 = tri[i2+1];
			vec3d n1 = (sm->vList[t1.v0].c-sm->vList[t1.v1].c)^(sm->vList[t1.v0].c-sm->vList[t1.v2].c);
			if (n1.squaredLength() > 0)	n1.normalize();
			vec3d n2 = (sm->vList[t2.v0].c-sm->vList[t2.v1].c)^(sm->vList[t2.v0].c-sm->vList[t2.v2].c);
			if (n2.squaredLength() > 0) n2.normalize();
			w += (1-(n1|n2));
			//n1.print();
			//n2.print();
		}
		if (w < best_weight) {
			best_weight = w;
			best_id = i1;
		}
		//std::cerr << i1 << " -> " << w << std::endl;
	}
	//std::cerr << "> " << best_id << " -> " << best_weight << std::endl;

	return best_id;
}

//******************************************************************************************
//*********************************** TRIANGLE *********************************************
//******************************************************************************************


Triangle::Triangle(void) {
	marker = false;
}

//******************************************************************************************

Triangle::Triangle(int v0, int v1, int v2, Edge *e0_, Edge *e1_, Edge *e2_, SimpleMesh* mesh)
{
	m_v0 = v0;
	m_v1 = v1;
	m_v2 = v2;
	e0 = e0_;
	e1 = e1_;
	e2 = e2_;
	m_mesh = mesh;
	marker = false;
}

//******************************************************************************************

Triangle::~Triangle(void)
{
}

//******************************************************************************************

#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
vec3f Triangle::getInterpolatedColor(const vec3d& barys) const {
	return m_mesh->vList[m_v0].color*(float)barys.x + m_mesh->vList[m_v1].color*(float)barys.y + m_mesh->vList[m_v2].color*(float)barys.z;
}
#endif


//******************************************************************************************

vec3d Triangle::getInterpolatedPos(const vec3d& barys) const {
	return m_mesh->vList[m_v0].c*(float)barys.x + m_mesh->vList[m_v1].c*(float)barys.y + m_mesh->vList[m_v2].c*(float)barys.z;
}

//******************************************************************************************

vec3d Triangle::getNormal() const {
	vec3d normal = (m_mesh->vList[m_v0].c - m_mesh->vList[m_v1].c)^(m_mesh->vList[m_v0].c - m_mesh->vList[m_v2].c);
	normal.normalize();
	return normal;
}

//******************************************************************************************

double Triangle::getArea() const {
	vec3d normal = (m_mesh->vList[m_v0].c - m_mesh->vList[m_v1].c)^(m_mesh->vList[m_v0].c - m_mesh->vList[m_v2].c);
	return normal.length()/2.0;
}

//******************************************************************************************

vec3d Triangle::barycentric(vec2d p) const {

	vec3d v0 = m_mesh->getVertex(m_v0);
	vec3d v1 = m_mesh->getVertex(m_v1);
	vec3d v2 = m_mesh->getVertex(m_v2);

	double A_tri = v1[0]*v2[1] + v2[0]*v0[1] + v0[0]*v1[1] - v0[0]*v2[1] - v2[0]*v1[1] - v1[0]*v0[1];
	double A1 =  p[0]*v2[1] + v2[0]*v0[1] + v0[0]*p[1] - v0[0]*v2[1] - v2[0]*p[1] - p[0]*v0[1];
	double A2 = v1[0]*p[1] + p[0]*v0[1] + v0[0]*v1[1] - v0[0]*p[1] - p[0]*v1[1] - v1[0]*v0[1];
	double A0 = v1[0]*v2[1] + v2[0]*p[1] + p[0]*v1[1] - p[0]*v2[1] - v2[0]*v1[1] - v1[0]*p[1];

	double u = A0/A_tri;
	double v = A1/A_tri;
	double w = A2/A_tri;

	if (fabs(A_tri) < 0.00000000001) {
		std::cerr << "Warning, very small triangle" << std::endl;
		std::cerr << "size: " << A_tri << std::endl;
		//exit(1);
	}

	return vec3d(u, v, w);
}

//******************************************************************************************

vec3d Triangle::barycentric_with_proj(vec3d& p) const {

	const vec3d& v0 = m_mesh->vList[m_v0].c;
	const vec3d& v1 = m_mesh->vList[m_v1].c;
	const vec3d& v2 = m_mesh->vList[m_v2].c;

	// project p onto triangle plane
	vec3d normal = (v0-v1)^(v0-v2);
	//double A_tri = normal.length();
	normal.normalize();

	p = p-normal*(normal|(p-v0));

	vec3d t2 = ((v0-v1)^(v0-p));
	double A2 = t2.length();
	vec3d t1 = ((v0-p)^(v0-v2));
	double A1 = t1.length();
	vec3d t0 = ((p-v1)^(p-v2));
	double A0 = t0.length();

	double sum = A0+A1+A2;

	if ((t2|normal) < 0) A2 *= -1;
	if ((t1|normal) < 0) A1 *= -1;
	if ((t0|normal) < 0) A0 *= -1;

	return vec3d(A0, A1, A2)/sum;
}

//******************************************************************************************

vec3d Triangle::barycentric(vec3d p) const {

	vec3d v0 = m_mesh->getVertex(m_v0);
	vec3d v1 = m_mesh->getVertex(m_v1);
	vec3d v2 = m_mesh->getVertex(m_v2);

	// project p onto triangle plane
	vec3d normal = (v0-v1)^(v0-v2);
	//double A_tri = normal.length();
	normal.normalize();

	p = p-normal*(normal|(p-v0));

	vec3d t2 = ((v0-v1)^(v0-p));
	double A2 = t2.length();
	vec3d t1 = ((v0-p)^(v0-v2));
	double A1 = t1.length();
	vec3d t0 = ((p-v1)^(p-v2));
	double A0 = t0.length();

	double sum = A0+A1+A2;

	if ((t2|normal) < 0) A2 *= -1;
	if ((t1|normal) < 0) A1 *= -1;
	if ((t0|normal) < 0) A0 *= -1;

	return vec3d(A0, A1, A2)/sum;
}

//******************************************************************************************

vec3d Triangle::getCenter() const {
	vec3d center = (m_mesh->vList[m_v0].c+m_mesh->vList[m_v1].c+m_mesh->vList[m_v2].c)/3;
	return center;
}

//******************************************************************************************

int& Triangle::operator[](unsigned int i) {
	if (i == 0)
		return m_v0;
	else if (i == 1)
		return m_v1;
	else if (i == 2)
		return m_v2;

	// Should NOT happen
	else
		exit(EXIT_FAILURE);
		return m_v0;
}

//******************************************************************************************

int Triangle::compare(Triangle other) {

	int matches = 0;
	for (int i1 = 0; i1 < 3; i1++)
		for (int i2 = 0; i2 < 3; i2++)
			if ((*this)[i1] == other[i2])
				matches++;

	return matches;
}

//******************************************************************************************

Triangle* Triangle::getTriangleNeighbor(int i) {
	Edge* e;

	if (i == 0) e = e0;
	else if (i == 1) e = e1;
	else e = e2;

	if (e->tList[0] != this) return e->tList[0];
	else {
		if (e->tList.size() == 1) return 0;
		else return e->tList[1];
	}
}

//******************************************************************************************

int Triangle::getOther(int v0, int v1) {
	if ((m_v0 != v0) && (m_v0 != v1))
		return m_v0;
	else if ((m_v1 != v0) && (m_v1 != v1))
		return m_v1;
	else
		return m_v2;
}

//******************************************************************************************

double Triangle::getAngleAtVertex(int v) {
	vec3d dir1, dir2;


	if (v == m_v0) {
		dir1 = m_mesh->vList[m_v0].c - m_mesh->vList[m_v1].c;
		dir2 = m_mesh->vList[m_v0].c - m_mesh->vList[m_v2].c;
	} else if (v == m_v1) {
		dir1 = m_mesh->vList[m_v1].c - m_mesh->vList[m_v0].c;
		dir2 = m_mesh->vList[m_v1].c - m_mesh->vList[m_v2].c;
	} else if (v == m_v2) {
		dir1 = m_mesh->vList[m_v2].c - m_mesh->vList[m_v0].c;
		dir2 = m_mesh->vList[m_v2].c - m_mesh->vList[m_v1].c;
	} else {
		return 0.0;
	}

	dir1.normalize();
	dir2.normalize();

	return acos(dir1|dir2);
}


//******************************************************************************************

void Triangle::print() {
	std::cerr << m_v0 << " " << m_v1 << " " << m_v2 << std::endl;
}


//******************************************************************************************

bool Triangle::isConsistent() const {
	std::set<int> vs;
	vs.insert(e0->v0);
	vs.insert(e0->v1);
	vs.insert(e1->v0);
	vs.insert(e1->v1);
	vs.insert(e2->v0);
	vs.insert(e2->v1);

	if (vs.size() != 3) {
		return false;
	}

	if (vs.find(m_v0) == vs.end())
		return false;
	if (vs.find(m_v1) == vs.end())
		return false;
	if (vs.find(m_v2) == vs.end())
		return false;

	if (e0->tList.size() == 2) if (e0->tList[0]==e0->tList[1])
		return false;
	if (e1->tList.size() == 2) if (e1->tList[0]==e1->tList[1])
		return false;
	if (e2->tList.size() == 2) if (e2->tList[0]==e2->tList[1])
		return false;

	return true;
}

//******************************************************************************************
//******************************** TRIANGLELIST ********************************************
//******************************************************************************************


//TriangleList::~TriangleList() {
//	for (unsigned int i1 = 0; i1 < size(); i1++) {
//		delete (*this)[i1];
//	}
//}

//******************************************************************************************

void TriangleList::insertTriangle(Triangle* t) {
	push_back(t);
}

//******************************************************************************************

bool TriangleList::removeTriangle(int v0, int v1, int v2) {

	Triangle t(v0, v1, v2, NULL, NULL, NULL, NULL);

	// Search the triangle:
	iterator tIter;

	for(tIter = begin(); tIter != end(); tIter++) {
		if(t.compare(*(*tIter)) == 3)
			break;
	}
	// Triangle is not contained in the list
	if(tIter == end())
		return false;
	else {
		erase(tIter);
		return true;
	}
}

//******************************************************************************************

UInt TriangleList::getTriangleIndex(int v0, int v1, int v2) {

	Triangle t(v0, v1, v2, NULL, NULL, NULL, NULL);

	// Search the triangle:
	iterator tIter;

	for (UInt i1 = 0; i1 < size(); i1++) {
		if (t.compare(*(*this)[i1]) == 3) return i1;
	}

	// Triangle is not contained in the list
	return 0;
}

//******************************************************************************************

Triangle* TriangleList::getTriangle(int v0, int v1, int v2) {

	Triangle t(v0, v1, v2, NULL, NULL, NULL, NULL);

	// Search the triangle:
	iterator tIter;

	for(tIter = begin(); tIter != end(); tIter++) {
		if(t.compare(*(*tIter)) == 3)
			break;
	}
	// Triangle is not contained in the list
	if(tIter == end())
		return NULL;
	else {
		return *tIter;
	}
}


//******************************************************************************************
//************************************* VERTEX *********************************************
//******************************************************************************************


Vertex::Vertex() {
#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
	param = vec2d(0,0);
#endif
	bool_flag = false;
}

//******************************************************************************************

Vertex::Vertex(vec3d c_) : c(c_) {
#ifndef USE_LIGHTWEIGHT_SIMPLE_MESH
	param = vec2d(0,0);
#endif
	bool_flag = false;
}

//******************************************************************************************

Vertex::~Vertex() {}


//******************************************************************************************
//************************************** EDGE **********************************************
//******************************************************************************************


Edge::Edge(){
	marker = false;
}

//******************************************************************************************

Edge::Edge(int v0_, int v1_, SimpleMesh* mesh) {
	m_mesh = mesh;
	if (v0_ < v1_) {
		v0 = v0_;
		v1 = v1_;
	} else {
		v1 = v0_;
		v0 = v1_;
	}

	marker = false;
	m_mesh = mesh;
}

//******************************************************************************************

Edge::~Edge(){}

//******************************************************************************************

double Edge::computeDihedralAngle() const {

	if (tList.size() != 2)
		return 0;

	vec3d n0 = tList[0]->getNormal();
	vec3d n1 = tList[1]->getNormal();

	//std::cerr << ((int)tList[0])%10000 << " " << ((int)tList[1])%10000 << " ";

	return 1-(n0|n1);
}

//******************************************************************************************

double Edge::computeDihedralAngles() const {

	if (tList.size() != 2)
		return 0;

	double sum = computeDihedralAngle();

	Triangle* t0 = tList[0];
	Triangle* t1 = tList[1];

	if (t0->e0 != this) sum += t0->e0->computeDihedralAngle();
	if (t0->e1 != this) sum += t0->e1->computeDihedralAngle();
	if (t0->e2 != this) sum += t0->e2->computeDihedralAngle();
	if (t1->e0 != this) sum += t1->e0->computeDihedralAngle();
	if (t1->e1 != this) sum += t1->e1->computeDihedralAngle();
	if (t1->e2 != this) sum += t1->e2->computeDihedralAngle();

	return sum;
}

//******************************************************************************************

void Edge::computeCotangentWeightsVERBOSE() const {

	if (tList.size() == 0) return;

	for (UInt i1 = 0; i1 < tList.size(); i1++) {

		int v2 = tList[i1]->getOther(v0, v1);
		vec3d vec0 = m_mesh->vList[v2].c-m_mesh->vList[v1].c;
		//vec0.normalize();
		vec3d vec1 = m_mesh->vList[v2].c-m_mesh->vList[v0].c;
		//vec1.normalize();

		//double cos_ = (vec0|vec1)/(vec0.length()*vec1.length());
		//double sin_ = (vec0^vec1).length()/(vec0.length()*vec1.length());

		std::cerr << (vec0|vec1)/(vec0^vec1).length() << std::endl;
	}

}

//******************************************************************************************

void Edge::computeCotangentWeights(bool uniform) {

	if (uniform) {
		cotangent_weight = 1;
		return;
	}
	
	cotangent_weight = 0;

	if (tList.size() == 0) return;

	for (UInt i1 = 0; i1 < tList.size(); i1++) {

		int v2 = tList[i1]->getOther(v0, v1);
		vec3d vec0 = m_mesh->vList[v2].c-m_mesh->vList[v1].c;
		//vec0.normalize();
		vec3d vec1 = m_mesh->vList[v2].c-m_mesh->vList[v0].c;
		//vec1.normalize();

		//double cos_ = (vec0|vec1)/(vec0.length()*vec1.length());
		//double sin_ = (vec0^vec1).length()/(vec0.length()*vec1.length());

		cotangent_weight += (vec0|vec1)/(vec0^vec1).length();
	}
	//cotangent_weight /= tList.size(); // Never uncomment this line unless you really knwo what you're doing!

}

//******************************************************************************************

void Edge::Swap() {

	//std::cerr << "v0: " << v0;
	//m_mesh->vList[v0].eList.print();
	//std::cerr << "v1: " << v1;
	//m_mesh->vList[v1].eList.print();

	//return;

	if (tList.size() != 2) return;

	Triangle* t0 = tList[0];
	vec3d n0 = t0->getNormal();
	Triangle* t1 = tList[1];
	vec3d n1 = t1->getNormal();
	int v2 = t0->getOther(v0, v1);
	int v3 = t1->getOther(v0, v1);

	// test if egde (v2, v3) exists:
	if (m_mesh->vList[v2].eList.getEdge(v2, v3) != NULL) {
		std::cerr << "Tetraeder" << std::endl;
		return;
	}
	if (t0->compare(*t1) == 3) {
		std::cerr << "Triangles are equal!!" << std::endl;
		return;
	}

	Edge* e1t0;		// Edge v0 v2
	e1t0 = m_mesh->vList[v0].eList.getEdge(v0, v2);
	Edge* e2t0;		// Edge v2 v1
	e2t0 = m_mesh->vList[v1].eList.getEdge(v2, v1);

	Edge* e1t1;		// Edge v0 v3
	e1t1 = m_mesh->vList[v0].eList.getEdge(v0, v3);
	Edge* e2t1;		// Edge v3 v1
	e2t1 = m_mesh->vList[v1].eList.getEdge(v3, v1);

	if (e1t0 == NULL) {
		std::cerr << "e1t0 is NULL" << std::endl;
		return;
	}
	if (e2t0 == NULL) {
		std::cerr << "e2t0 is NULL" << std::endl;
		return;
	}
	if (e1t1 == NULL) {
		std::cerr << "e1t1 is NULL" << std::endl;
		return;
	}
	if (e2t1 == NULL) {
		std::cerr << "e2t1 is NULL" << std::endl;
		return;
	}

	vec3d n2;
	if (e1t0->tList.size() != 2) {
		//std::cerr << "Warning: e1t0" << std::endl;
		n2 = vec3d(0.0);
	} else {
		n2 = (e1t0->tList[0] != t0) ? e1t0->tList[0]->getNormal() : e1t0->tList[1]->getNormal();
	}
	vec3d n3;
	if (e2t0->tList.size() != 2) {
		//std::cerr << "Warning: e2t0" << std::endl;
		n3 = vec3d(0.0);
	} else {
		n3 = (e2t0->tList[0] != t0) ? e2t0->tList[0]->getNormal() : e2t0->tList[1]->getNormal();
	}
	vec3d n4;
	if (e1t1->tList.size() != 2) {
		//std::cerr << "Warning: e1t1" << std::endl;
		n4 = vec3d(0.0);
	} else {
		n4 = (e1t1->tList[0] != t1) ? e1t1->tList[0]->getNormal() : e1t1->tList[1]->getNormal();
	}
	vec3d n5;
	if (e2t1->tList.size() != 2) {
		//std::cerr << "Warning: e2t1" << std::endl;
		n5 = vec3d(0.0);
	} else {
		n5 = (e2t1->tList[0] != t1) ? e2t1->tList[0]->getNormal() : e2t1->tList[1]->getNormal();
	}

	double current_curv = 1-(n0|n1) + 1-(n0|n2) + 1-(n0|n3) + 1-(n1|n4) + 1-(n1|n5);

	vec3d n0_new = (m_mesh->vList[v0].c-m_mesh->vList[v3].c)^(m_mesh->vList[v0].c-m_mesh->vList[v2].c);
	vec3d n1_new = (m_mesh->vList[v1].c-m_mesh->vList[v2].c)^(m_mesh->vList[v1].c-m_mesh->vList[v3].c);
	n0_new.normalize();
	n1_new.normalize();

	//if ((n0_new|n0) < 0) {
	//	return;
	//	n0_new *= -1;
	//	std::cerr << "oops n0" << std::endl;
	//}
	//if ((n1_new|n1) < 0) {
	//	return;
	//	n1_new *= -1;
	//	std::cerr << "oops n1" << std::endl;
	//}

	double curv_after_swap = 1-(n0_new|n1_new) + 1-(n0_new|n2) + 1-(n0_new|n4) + 1-(n1_new|n3) + 1-(n1_new|n5);

	if ((curv_after_swap - current_curv) > 0.000001) // don't swap.
		return;

	std::cerr << "Curv new:	\t" << curv_after_swap << " Curv old \t" << current_curv << " \t" << computeDihedralAngles() << std::endl;

	// Update triangles
	t0->set_v0(v0); t0->set_v1(v3); t0->set_v2(v2);
	t1->set_v0(v1); t1->set_v1(v2); t1->set_v2(v3);

	// Update triangle-lists
	if (e2t0->tList[0] == t0) e2t0->tList[0] = t1;
	if (e2t0->tList[1] == t0) e2t0->tList[1] = t1;
	if (e1t1->tList[0] == t1) e1t1->tList[0] = t0;
	if (e1t1->tList[1] == t1) e1t1->tList[1] = t0;

	//if (e2t0->tList[0]->compare(*(e2t0->tList[1])) == 3) {
	//	std::cerr << "OOPS: e2t0Triangles are equal!!" << std::endl;
	//	e2t0->tList.print();
	//	return;
	//}
	//if (e1t1->tList[0]->compare(*(e1t1->tList[1])) == 3) {
	//	std::cerr << "OOPS: e1t1Triangles are equal!!" << std::endl;
	//	e1t1->tList.print();
	//	return;
	//}

	// Update edges
	int old_v0 = v0; int old_v1 = v1;
	this->v0 = v2; this->v1 = v3;
	//this->print();
	// Update edges in vLists:
	m_mesh->vList[old_v0].eList.removeEdge(this);
	m_mesh->vList[old_v1].eList.removeEdge(this);
	m_mesh->vList[v2].eList.insertEdge(this);
	m_mesh->vList[v3].eList.insertEdge(this);
	// Update edges in triangles
	t0->e0 = e1t0; t0->e1 = this; t0->e2 = e1t1;
	t1->e0 = e2t0; t1->e1 = e2t1; t1->e2 = this;

	std::cerr << "after: \t" << computeDihedralAngles() << std::endl;

	if (!t0->isConsistent() || !t1->isConsistent()) {
		std::cerr << "Inconsistent triangle" << std::endl;
	}

	// Tell other dudes to swap.
	//std::cerr << "swapped" << std::endl;
	e1t0->Swap();
	e2t0->Swap();
	e1t1->Swap();
	e2t1->Swap();
}

//******************************************************************************************

bool Edge::operator==(const Edge& other) {
	if (v0 == other.v0 && v1 == other.v1)
		return true;
	else
		return false;
}

//******************************************************************************************

double Edge::getLength() const {
	return (m_mesh->vList[v0].c).dist(m_mesh->vList[v1].c);
}

//******************************************************************************************

void Edge::print() {
	std::cerr << "Edge:" << std::endl;
	std::cerr << v0 << " " << v1 << std::endl;
	std::cerr << "Adjacent triangles:" << std::endl;
	for (UInt i1 = 0; i1 < (UInt)tList.size(); i1++) {
		tList[i1]->print();
	}
}

//******************************************************************************************

UInt Edge::getIndexOfTriangle(UInt i) {

	Triangle* t0 = tList[i];
	UInt pos = m_mesh->tList.getTriangleIndex(t0->getV(0), t0->getV(1), t0->getV(2));

	return pos;
}

//******************************************************************************************
//********************************** EDGE LIST *********************************************
//******************************************************************************************


Edge* EdgeList::getEdge(int v0, int v1) {

	if (size() == 0)
		return NULL;

	iterator it;

	for (it = begin(); it < end(); it++)
		if (v0 == (*it)->v0 && v1 == (*it)->v1 || v1 == (*it)->v0 && v0 == (*it)->v1)
			return (*it);

	return NULL;
}

//******************************************************************************************

void EdgeList::print() {

	iterator it;
	std::cerr << "EdgeList:" << std::endl;
	for (it = begin(); it < end(); it++)
		(*it)->print();
}

//******************************************************************************************

void EdgeList::insertEdge(Edge* e) {
	push_back(e);
}

//******************************************************************************************

bool EdgeList::removeEdge(Edge* e) {

	if (size() == 0)
		return false;

	return removeEdge(e->v0, e->v1);
}

//******************************************************************************************

bool EdgeList::removeEdge(int v0, int v1) {

	if (size() == 0)
		return false;

	iterator it;

	for (it = begin(); it < end(); it++)
		if (v0 == (*it)->v0 && v1 == (*it)->v1 || v1 == (*it)->v0 && v0 == (*it)->v1)
			break;

	// Edge is not contained in the list
	if(it == end()) {
		std::cerr << "EdgeList::removeEdge(...): Edge not contained" << std::endl;
		return false;
	} else {
		// std::cerr << "Edge " << v0 << ", " << v1 << " removed" << std::endl;
 		erase(it);
		return true;
	}
}


//******************************************************************************************
//********************************** SIMPLE MESH *******************************************
//******************************************************************************************

SimpleMesh::SimpleMesh() {
	m_numV = 0;
	m_numT = 0;
	boundary_vertices_computed = false;
	hasTextureCoordinates = false;
	hasColors = false;
}

//******************************************************************************************

SimpleMesh::~SimpleMesh(void) {


	//std::cerr << "SimpleMesh destructor called...";

	// Delete all Triangles
	for (unsigned int i1 = 0; i1 < tList.size(); i1++)
		delete tList[i1];
	tList.clear();

	vList.clear();

	//std::cerr << "Deleting edges" << std::endl;
	// Delete all Edges
	for (unsigned int i1 = 0; i1 < eList.size(); i1++)
		delete eList[i1];
	eList.clear();
	//std::cerr << "Deleting edges done" << std::endl;

	//std::cerr << "destroyed" << std::endl;

}


//******************************************************************************************

std::vector<int> SimpleMesh::getBoundary() {
	std::vector<int> boundary;

	//  Search for first edge with only one adjacent triangle:
	Edge* e = NULL;

	for (UInt i1 = 0; i1 < eList.size(); i1++) eList[i1]->marker = false;

	int start = 0;

	do {

		//std::cerr << "Start at: " << start << std::endl;

		for (unsigned int i1 = start; i1 < eList.size(); i1++) {
			e = eList[i1];
			start = i1;
			if (e->tList.size() == 1)
				break;
		}
		start++;

		if (start >= (int)eList.size()) break;

		bool doesntmatter;

		boundary = getBoundaryStartAtEdge(e, doesntmatter);
		if (boundary.size() > 0)
			return boundary;

	} while (true);

	
	return boundary;

}

//******************************************************************************************

std::vector<int> SimpleMesh::getBoundaryStartAtEdge(Edge* e, bool& is_included) {

	//Edge* org_edge = e;
	e->marker = true;

	std::vector<int> boundary;
	std::set<int> boundary_set;

	int start_v = e->v0;
	int next_v = e->v1;
	Triangle* t = e->tList[0];
	if ((t->v0() == e->v1 && t->v1() == e->v0) || (t->v1() == e->v1 && t->v2() == e->v0) || (t->v2() == e->v1 && t->v0() == e->v0)) {
		start_v = e->v1;
		next_v = e->v0;
	}



	boundary.push_back(start_v);
	boundary.push_back(next_v);
	boundary_set.insert(start_v);
	boundary_set.insert(next_v);

	while(next_v != start_v) {
		e = NULL;
		unsigned int i1;
		EdgeList current_eList = vList[next_v].eList;
		for (i1 = 0; i1 < current_eList.size(); i1++) {
			if (current_eList[i1]->marker == false && current_eList[i1]->tList.size() == 1 && current_eList[i1]->v0 != boundary[boundary.size()-2] && current_eList[i1]->v1 != boundary[boundary.size()-2]) {
				e = current_eList[i1];
				break;
			}
		}

		if (e == NULL) {
			std::cerr << "Ungeschlossener Randzyklus gefunden!" << std::endl;
			std::cerr << "Boundary: ";
			for (unsigned int i1 = 0; i1 < boundary.size(); i1++)
				std::cerr << boundary[i1] << " ";
			std::cerr << std::endl;
			return std::vector<int>();
		}

		e->marker = true;

		next_v = (current_eList[i1]->v0 == boundary[boundary.size()-1]) ? e->v1 : e->v0;
		boundary.push_back(next_v);

		if (boundary_set.find(next_v) != boundary_set.end()) {
			break;
		}
		boundary_set.insert(next_v);
	}

	if (boundary[0] == boundary[boundary.size()-1]) {
		is_included = true;
		return boundary;
	}

	//std::cerr << "Boundary: ";
	//for (unsigned int i1 = 0; i1 < boundary.size(); i1++)
	//	std::cerr << boundary[i1] << " ";
	//std::cerr << std::endl;

	// We have a cycle:
	std::vector<int> cycle;
	cycle.push_back(next_v);
	UInt i1_save = 0;
	for (UInt i1 = (UInt)boundary.size()-2; i1 > 0; i1--) {
		i1_save = i1;

		int current = boundary[i1];
		cycle.push_back(current);
		if (current == next_v)
			break;
	}

	//std::cerr << "cycle: ";
	//for (unsigned int i1 = 0; i1 < cycle.size(); i1++)
	//	std::cerr << cycle[i1] << " ";
	//std::cerr << std::endl;

	//std::cerr << "i1_save: " << i1_save << std::endl;

	// Unmark remaining:
	//std::cerr << "Unmark: " << boundary[i1_save];
	for (UInt i1 = i1_save; i1 > 0; i1--) {
		int current = boundary[i1];
		int last = boundary[i1-1];

		Edge* e = vList[current].eList.getEdge(current, last);
		e->marker = false;

		//std::cerr << last << " ";
	}

	//exit(1);

	//std::cerr << "Cycle found" << std::endl;
	//std::cerr << "Boundary: ";
	//for (unsigned int i1 = 0; i1 < boundary.size(); i1++)
	//	std::cerr << boundary[i1] << " ";
	//std::cerr << std::endl;
	//std::cerr << "Cycle: ";
	//for (unsigned int i1 = 0; i1 < cycle.size(); i1++)
	//	std::cerr << cycle[i1] << " ";
	//std::cerr << std::endl;

	is_included = false;
	return cycle;
}

//******************************************************************************************

Edge* SimpleMesh::insertEdge(int v0, int v1) {

	if (v0 < 0 || v0 >= (int)vList.size() || v1 < 0 || v1 >= (int)vList.size())
		std::cerr << "Fehler bei SimpleMesh::insertEdge(" << v0 << ", " << v1 << ")! Vertex index out of range!" << std::endl;

	Edge* e = vList[v0].eList.getEdge(v0, v1);
	if (e != NULL) {
		//e->print();
		//std::cerr << "Fehler bei SimpleMesh::insertEdge(" << v0 << ", " << v1 << "): Kante existiert bereits" << std::endl;
		return e;
	}

	e = new Edge(v0, v1, this);

	// Insert edge in verticees and global list
	vList[v0].eList.insertEdge(e);
	vList[v1].eList.insertEdge(e);
	eList.insertEdge(e);

	return e;
}

//******************************************************************************************

Triangle* SimpleMesh::insertTriangle(int v0, int v1, int v2) {

	//std::cerr << "SimpleMesh::insertTriangle(" << v0 << ", " << v1 << ", " << v2 << ")" << std::endl;

	if (v0 < 0 || v0 >= (int)vList.size() || v1 < 0 || v1 >= (int)vList.size() || v2 < 0 || v2 >= (int)vList.size()) {
#ifndef SIMPLEMESH_SHUT_UP
		std::cerr << "Fehler bei SimpleMesh::insertTriangle(" << v0 << ", " << v1 << ", " << v2 << ")! Vertex index out of range!" << std::endl;
#endif
	}

	if ((v0 == v1) || (v0 == v2) || (v1 == v2)) {
#ifndef SIMPLEMESH_SHUT_UP
		std::cerr << "Fehler bei SimpleMesh::insertTriangle(" << v0 << ", " << v1 << ", " << v2 << ")! Inserted degenerated triangle!" << std::endl;
#endif
		return 0;
	}

	Edge* e0 = vList[v1].eList.getEdge(v1, v2);
	if(e0 == NULL)
		e0 = insertEdge(v1, v2);

	Edge* e1 = vList[v0].eList.getEdge(v2, v0);
	if(e1 == NULL)
		e1 = insertEdge(v2, v0);

	Edge* e2 = vList[v1].eList.getEdge(v0, v1);
	if(e2 == NULL)
		e2 = insertEdge(v0, v1);


	if ((e0->tList.size() > 2) || (e1->tList.size() > 2) || (e2->tList.size() > 2)) {
#ifndef SIMPLEMESH_SHUT_UP
		std::cerr << "Fehler bei SimpleMesh::insertTriangle(" << v0 << ", " << v1 << ", " << v2 << ")! Inserted triangle at full edge!" << std::endl;
#endif
		return 0;
	}

	Triangle* t = new Triangle(v0, v1, v2, e0, e1, e2, this);
	tList.insertTriangle(t);
	e0->tList.insertTriangle(t);
	e1->tList.insertTriangle(t);
	e2->tList.insertTriangle(t);



	m_numT++;


	return t;
}

//******************************************************************************************

void SimpleMesh::insertVertex(double x, double y, double z) {
	vList.push_back(Vertex(vec3d(x,y,z)));
	m_numV++;

	// Calculate extend of the Mesh
	if(m_numV == 1) {
		min[0] = x;
		max[0] = x;
		min[1] = y;
		max[1] = y;
		min[2] = z;
		max[2] = z;
	} else {
		if(x < min[0])
			min[0] = x;
		if(x > max[0])
			max[0] = x;
		if(y < min[1])
			min[1] = y;
		if(y > max[1])
			max[1] = y;
		if(z < min[2])
			min[2] = z;
		if(z > max[2])
			max[2] = z;
	}
}

//******************************************************************************************

void SimpleMesh::print() {

	std::cerr << "====================================================" << std::endl;
	std::cerr << "========            SIMPLE MESH           ==========" << std::endl;
	std::cerr << "====================================================" << std::endl;
	std::cerr << "Points" << std::endl;
	for (int i1 = 0; i1 < m_numV; i1++) {
		// Print coordinates
		vList[i1].c.print();
	}
	std::cerr << "Triangles" << std::endl;
	for (int i1 = 0; i1 < m_numT; i1++) {
		// Print triangles
		tList[i1]->print();
	}
	std::cerr << "Edges" << std::endl;
	for (int i1 = 0; i1 < (int)eList.size(); i1++) {
		// Print triangles
		std::cerr << eList[i1]->v0 << " " << eList[i1]->v1 << " numTris: " << (int)eList[i1]->tList.size() << std::endl;
	}
}

//******************************************************************************************

void SimpleMesh::computeBoundaryVertices() {

	for (unsigned int i1 = 0; i1 < vList.size(); i1++) vList[i1].is_boundary_point = false;


	for (unsigned int i1 = 0; i1 < eList.size(); i1++) {
		Edge* e = eList[i1];
		if (e->tList.size() < 2) { // It's a boundary edge
			vList[e->v0].is_boundary_point = true;
			vList[e->v1].is_boundary_point = true;
		}
	}

	boundary_vertices_computed = true;
}

//******************************************************************************************

void SimpleMesh::removeTriangle(Triangle* t) {
	int v0 = t->v0();
	int v1 = t->v1();
	int v2 = t->v2();

	Edge* e0 = t->e0;
	Edge* e1 = t->e1;
	Edge* e2 = t->e2;

	if (!e0->tList.removeTriangle(v0, v1, v2) || !e1->tList.removeTriangle(v0, v1, v2) || !e2->tList.removeTriangle(v0, v1, v2))
		std::cerr << "SimpleMesh::removeTriangle(...): Triangle could not be removed from edges!)" << std::endl;

	if (!tList.removeTriangle(v0, v1, v2))
		std::cerr << "SimpleMesh::removeTriangle(...): Triangle could not be removed from list!)" << std::endl;

	delete t;

	if (e0->tList.size() == 0) removeEdge(e0);
	if (e1->tList.size() == 0) removeEdge(e1);
	if (e2->tList.size() == 0) removeEdge(e2);

	m_numT--;
}

//******************************************************************************************

void SimpleMesh::removeTriangleNoPhysical(Triangle* t) {
	int v0 = t->v0();
	int v1 = t->v1();
	int v2 = t->v2();

	Edge* e0 = t->e0;
	Edge* e1 = t->e1;
	Edge* e2 = t->e2;

	if (!e0->tList.removeTriangle(v0, v1, v2) || !e1->tList.removeTriangle(v0, v1, v2) || !e2->tList.removeTriangle(v0, v1, v2))
		std::cerr << "SimpleMesh::removeTriangle(...): Triangle could not be removed from edges!)" << std::endl;

	t->marker = true;

	if (e0->tList.size() == 0) removeEdge(e0);
	if (e1->tList.size() == 0) removeEdge(e1);
	if (e2->tList.size() == 0) removeEdge(e2);
}

//******************************************************************************************

void SimpleMesh::cleanUpMarkedTriangles() {
	for (std::vector<Triangle*>::iterator it = tList.begin(); it != tList.end(); ) {
		Triangle* t = (*it);
		if (t->marker) {
			it = tList.erase(it);
			delete t;
		} else {
			it++;
		}
	}
	m_numT = (int)tList.size();
}

//******************************************************************************************

void SimpleMesh::removeVertex(int id) {
	std::set<Triangle*> trisAtVertex;

	std::vector<int> fan;
	for (UInt i1 = 0; i1 < vList[id].eList.size(); i1++) {
		Edge* e = vList[id].eList[i1];
		for (UInt i2 = 0; i2 < e->tList.size(); i2++) {
			Triangle* t = e->tList[i2];
			trisAtVertex.insert(t);
		}
		fan.push_back(e->getOther(id));
	}

	for (std::set<Triangle*>::iterator it = trisAtVertex.begin(); it != trisAtVertex.end(); it++) {
		//removeTriangleNoPhysical(*it);
		removeTriangle(*it);
	}
	vList[id].eList.clear();

	// Edges are automatically removed during 'removeTriangle'.

	//removeVertexPhysically(id);
}

//******************************************************************************************

void SimpleMesh::removeVertexPhysically(int id) {

	for (std::vector<Edge*>::iterator it = eList.begin(); it != eList.end(); it++) {
		Edge* e = (*it);
		if ((e->v0 == id) || (e->v0 == id)) {
			std::cerr << "Error in removeVertexPhysically: reference to vertex " << id << " left in edge  " << e->v0 << "-" << e->v1 << " !" <<  std::endl;
			std::cerr << e->tList.size() << std::endl;
		}

		if (e->v0 > id) e->v0 -= 1;
		if (e->v1 > id) e->v1 -= 1;
	}

	for (UInt i1 = 0; i1 < tList.size(); i1++) {
		Triangle* t = tList[i1];
		if ((t->m_v0 == id) || (t->m_v1 == id) || (t->m_v2 == id)) {
			std::cerr << "Error in removeVertexPhysically: reference to vertex " << id << " left in triangle  " << t->m_v0 << "-" << t->m_v1 << "-" << t->m_v2 << " !" << std::endl;
		}

		if (t->m_v0 > id) t->m_v0 -= 1;
		if (t->m_v1 > id) t->m_v1 -= 1;
		if (t->m_v2 > id) t->m_v2 -= 1;
	}

	vList.erase(vList.begin() + id);
	m_numV--;

}

//******************************************************************************************

void SimpleMesh::removeEdge(int v0, int v1) {

	Edge* e = vList[v0].eList.getEdge(v0, v1);
	removeEdge(e);

}

//******************************************************************************************

void SimpleMesh::removeEdge(Edge* e) {

	if (e == NULL) return;

	// Remove all triangles.
	std::vector<Triangle*> tris;
	for (UInt i1 = 0; i1 < e->tList.size(); i1++) {
		tris.push_back(e->tList[i1]);
	}
	for (UInt i1 = 0; i1 < tris.size(); i1++) {
		removeTriangle(tris[i1]);
	}

	const int& v0 = e->v0;
	const int& v1 = e->v1;

	// Remove edge from vertices
	vList[v0].eList.removeEdge(v0,v1);
	vList[v1].eList.removeEdge(v0,v1);
	// Remove edge from global list
	eList.removeEdge(v0,v1);
}

//******************************************************************************************

void SimpleMesh::cleanValence0Vertices(std::vector<edgeCollapsInfo>* eColInfo) {
	int* newIndices = new int[m_numV];

	int validVertCnt = 0;
	int vertexID = 0;
	for (std::vector<Vertex>::iterator it = vList.begin(); it != vList.end(); ) {
		newIndices[vertexID] = validVertCnt;
		vertexID++;
		if (it->getNumN() > 0) {
			validVertCnt++;
			it++;
		} else  {
			it = vList.erase(it);
		}
	}

	// Now update all triangles
	for (UInt i1 = 0; i1 < (UInt)m_numT; i1++) {
		Triangle* t = tList[i1];
		t->m_v0 = newIndices[t->m_v0];
		t->m_v1 = newIndices[t->m_v1];
		t->m_v2 = newIndices[t->m_v2];
	}

	// Now update all edges
	for (UInt i1 = 0; i1 < eList.size(); i1++) {
		Edge* e = eList[i1];
		e->v0 = newIndices[e->v0];
		e->v1 = newIndices[e->v1];
	}

	m_numV = (int)vList.size();


	if (eColInfo != NULL) {
		for (UInt i1 = 0; i1 < eColInfo->size(); i1++) {
			(*eColInfo)[i1].tri_corner_0 = newIndices[(*eColInfo)[i1].tri_corner_0];
			for (UInt i2 = 0; i2 < (*eColInfo)[i1].fan_edges.size(); i2++) {
				(*eColInfo)[i1].fan_edges[i2].x = newIndices[(*eColInfo)[i1].fan_edges[i2].x];
				(*eColInfo)[i1].fan_edges[i2].y = newIndices[(*eColInfo)[i1].fan_edges[i2].y];
			}
			(*eColInfo)[i1].basis.x = newIndices[(*eColInfo)[i1].basis.x];
			(*eColInfo)[i1].basis.y = newIndices[(*eColInfo)[i1].basis.y];
			(*eColInfo)[i1].basis.z = newIndices[(*eColInfo)[i1].basis.z];
		}
	}

	delete[] newIndices;
}

//******************************************************************************************


void SimpleMesh::computeOBB() {
	vec3d center;

	for (int i1 = 0; i1 < m_numV; i1++) {
		center += vList[i1].c;
	}
	center /= m_numV;

	SquareMatrixND<vec3d> CV;
	// Build covariance Matrix:
	for(int i1 = 0; i1 < m_numV; i1++) {
		vec3d diff = vList[i1].c - center;
		for(int i2 = 0; i2 < 3; i2++) {
			for(int i3 = i2; i3 < 3; i3++) {
				CV(i2, i3) += (diff[i2] * diff[i3]);
			}
		}
	}
	CV.symmetrize();
	SquareMatrixND<vec3d> ESystem;
	vec3d evals = CV.calcEValuesAndVectorsCORRECT(ESystem);
	// sort:

	vec3d base_x = ESystem.getColI(0);
	vec3d base_y = ESystem.getColI(1);
	vec3d base_z = ESystem.getColI(2);

	// Compute new center:
	vec3d max(-1000000, -1000000, -1000000);
	vec3d min(1000000, 1000000, 1000000);
	for(int i1 = 0; i1 < m_numV; i1++) {
		vec3d transp = center - vList[i1].c;
		double x = base_x|transp;
		double y = base_y|transp;
		double z = base_z|transp;

		if (x < min[0])
			min[0] = x;
		if (x > max[0])
			max[0] = x;
		if (y < min[1])
			min[1] = y;
		if (y > max[1])
			max[1] = y;
		if (z < min[2])
			min[2] = z;
		if (z > max[2])
			max[2] = z;
	}

	boundingbox.center = center;
	boundingbox.center -= base_x*(min[0]+max[0])/2;
	boundingbox.center -= base_y*(min[1]+max[1])/2;
	boundingbox.center -= base_z*(min[2]+max[2])/2;

    boundingbox.dx = (max[0]-min[0])/2;
    boundingbox.dy = (max[1]-min[1])/2;
    boundingbox.dz = (max[2]-min[2])/2;
	boundingbox.base_x = base_x;
	boundingbox.base_y = base_y;
	boundingbox.base_z = base_z;

	// Store BoundingBox basis:
	boundingbox.basis[0] = base_x[0];
	boundingbox.basis[1] = base_x[1];
	boundingbox.basis[2] = base_x[2];
	boundingbox.basis[3] = base_y[0];
	boundingbox.basis[4] = base_y[1];
	boundingbox.basis[5] = base_y[2];
	boundingbox.basis[6] = base_z[0];
	boundingbox.basis[7] = base_z[1];
	boundingbox.basis[8] = base_z[2];
}


//******************************************************************************************

vec3d SimpleMesh::getFanCentroid(int v) {

	if (vList[v].eList.size() == 0)	//isolated vertex
		return vec3d(0,0,0);

	vec3d centroid;

	for (unsigned int i1 = 0; i1 < vList[v].eList.size(); i1++) {
		centroid += vList[vList[v].eList[i1]->getOther(v)].c;
	}

	centroid /= (double)vList[v].eList.size();

	return centroid;
}

//******************************************************************************************
//******************************************************************************************
//******************************************************************************************


void UmbrellaSmooth(SimpleMesh* sm) {
	vec3d* new_verts = new vec3d[sm->getNumV()];
	sm->computeBoundaryVertices();
	for (UInt i1 = 0; i1 < (UInt)sm->getNumV(); i1++) {
		Vertex v = sm->vList[i1];
		if (v.is_boundary_point) {
			new_verts[i1] = v.c;
		} else { // Umbrella
			for (UInt i2 = 0; i2 < v.getNumN(); i2++) {
				UInt neighbor = v.getNeighbor(i2, i1);
				new_verts[i1] += sm->vList[neighbor].c;
			}
			new_verts[i1] /= v.getNumN();
		}
	}
	for (UInt i1 = 0; i1 < (UInt)sm->getNumV(); i1++) {
		sm->vList[i1].c = new_verts[i1];
	}
	delete[] new_verts;
}
