#pragma once

#include "util.h"

#include <SCAPE/SCAPE.h>

//
// conversions: Vec3d -> GLM
//
inline glm::vec3 toGLM(const vec3d& I) {
	return glm::vec3(I.x, I.y, I.z);
}
inline glm::vec3 toGLM(const Vertex& I) {
	return glm::vec3(I.c.x, I.c.y, I.c.z);
}
inline glm::mat3 toGLM(const SquareMatrixND<vec3d>& I) {
	glm::mat3 O = glm::identity<glm::mat3>();
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			O[x][y] = I(x, y);
		}
	}
	return O;
}
inline std::vector<glm::vec3> toGLM(const std::vector<vec3d>& I) {
	std::vector<glm::vec3> O;
	for (auto& i : I) {
		O.push_back(toGLM(i));
	}
	return O;
}
inline std::vector<glm::vec3> toGLM(const SimpleMesh& I) {
	std::vector<glm::vec3> O;
	for (auto& i : I.vList) {
		O.push_back(toGLM(i));
	}
	return O;
}

//
// conversions: GLM -> Vec3d
//
inline vec3d toVec3d(const glm::vec3& I) {
	return vec3d(I.x, I.y, I.z);
}
inline std::vector<vec3d> toVec3d(const std::vector<glm::vec3>& I) {
	std::vector<vec3d> O;
	for (auto& i : I) {
		O.push_back(toVec3d(i));
	}
	return O;
}


//
// SimpleMesh extension
//
inline glm::vec3 getCenterOfGravity(SimpleMesh * sm) {
	sm->computeBoundaryVertices();
	vec3d c = (sm->max - sm->min) / 2.0 + sm->min;

	return glm::vec3(c.x, c.y, c.z);
}
inline void transformMesh(SimpleMesh* sm, const glm::mat4& M) {
	for (int vertexID = 0; vertexID < sm->getNumV(); vertexID++) {
		glm::vec4 c = M * glm::vec4(toGLM(sm->vList.at(vertexID).c), 1.0);
		sm->vList.at(vertexID).c = toVec3d(glm::vec3(c));
	}
}
inline std::vector<glm::vec3> computeVertexNormals(SimpleMesh* sm) {
	std::vector<glm::vec3> normals(sm->vList.size(), glm::vec3(0.0));
	for (int i = 0; i < sm->getNumT(); i++) {
		Triangle* t = sm->tList.at(i);
		vec3d n = t->getNormal();
		glm::vec3 normal(n.x, n.y, n.z);

		normals.at(t->v0()) = normal;
		normals.at(t->v1()) = normal;
		normals.at(t->v2()) = normal;
	}

	return normals;
}
inline void transformMesh(SimpleMesh* sm, RigidTransformWithScale<vec3d> trans) {
	for (int vertexID = 0; vertexID < sm->vList.size(); vertexID++) {
		sm->vList.at(vertexID).c = trans.vecTrans(sm->vList.at(vertexID).c);
	}
}


//
// SCAPE
//
class SCAPE {
private:
	bool init = false;

	SCAPE() {};
	SCAPE(const SCAPE&) = delete;
	SCAPE& operator=(const SCAPE&) = delete;
public:
	static SCAPE& GetInstance() {
		static SCAPE instance;
		return instance;
	}

	SCAPE_Model scapeModel;

	inline int NumPoseParams() {
		return numParts * 3;
	}
	inline int NumShapeParams() {
		return scapeModel.pcabasis.numVecs;
	}
	inline int NumParams() {
		return NumPoseParams() + NumShapeParams();
	}

	void InitSCAPEModel(const std::string& binFile, const std::string& dataDirectory);
	void CheckInit();
	SimpleMesh CreateMeshFromSCAPEParams(const std::vector<double>& scapeParameters);
	std::vector<vec3d> ProjectToScapeSpace(const std::vector<glm::vec3>& vertices);
	std::vector<glm::vec3> FitScapeModel(const std::vector<glm::vec3>& vertices);
};
