#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <iterator>

#include <render/global.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/intersect.hpp>


//
// my data container structs
//
namespace Data {
	struct Mesh {
		std::vector<glm::vec3> vertices;
		std::vector<glm::vec3> normals;
		std::vector<uint> indices;
	};

	struct Marker {
		std::string name;
		glm::vec3 globalPosition;
		std::string attachedBody;
		glm::vec3 localPosition;
	};

	struct Body {
		std::string name;
		glm::mat4 globalTransform;
		glm::mat4 localTransform; // is always the identity
	};

	struct Joint {
		std::string name;
		std::string childName;
		std::string parentName;
		glm::mat4 globalTransformChild;
		glm::mat4 globalTransformParent;
		glm::mat4 localTransform; // is always the identity
	};

	struct Muscle {
		std::string name;
		std::vector<glm::vec3> lineSet;
		double activation;
		double maxForce;
	};
};

//
// usefull functions
//
inline std::string to_zero_lead(const int number, const int precision) {
	std::ostringstream oss;
	oss << std::setw(precision) << std::setfill('0') << number;
	return oss.str();
}

template<typename K, typename V>
inline bool FindKeyInMap(const K& key, const std::map<K, V>& map) {
	return map.find(key) != map.end();
}



// convert Data:: -> PointVector
inline std::vector<glm::vec3> toPointVector(const std::map<std::string, Data::Marker>& I) {
	std::vector<glm::vec3> O;
	for (const auto& i : I) {
		O.push_back(i.second.globalPosition);
	}
	return O;
}
inline std::vector<glm::vec3> toPointVector(const std::map<std::string, Data::Body>& I) {
	std::vector<glm::vec3> O;
	for (const auto& i : I) {
		O.push_back(glm::vec3(i.second.globalTransform[3]));
	}
	return O;
}

// convert Data:: -> PointMap
inline std::map<std::string, glm::vec3> toPointMap(const std::map<std::string, Data::Marker>& I) {
	std::map<std::string, glm::vec3> O;
	for (const auto& i : I) {
		O[i.first] = i.second.globalPosition;
	}
	return O;
}
inline std::map<std::string, glm::vec3> toPointMap(const std::map<std::string, Data::Body>& I) {
	std::map<std::string, glm::vec3> O;
	for (const auto& i : I) {
		O[i.first] = glm::vec3(i.second.globalTransform[3]);
	}
	return O;
}



//
// std::cout ostream <<
//
// map
template<typename K, typename V>
inline std::ostream& operator<<(std::ostream& out, const std::map<K, V>& m) {
	for (const std::pair<K, V>& p : m) {
		out << "{" << p.first << ": " << p.second << "}\n";
	}
	return out;
}
// vector
template <typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}

// GLM
inline std::ostream& operator<< (std::ostream& out, const glm::vec2& v) {
	return out << "(" << v.x << ", " << v.y << ")";
}
inline std::ostream& operator<< (std::ostream& out, const glm::vec3& v) {
	return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}
inline std::ostream& operator<< (std::ostream& out, const glm::vec4& v) {
	return out << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
}
inline std::ostream& operator<< (std::ostream& out, const glm::mat3& m) {
	out << std::endl;
	out << "(" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << "," << std::endl;
	out << " " << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << "," << std::endl;
	out << " " << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << ")";
	return out;
}
inline std::ostream& operator<< (std::ostream& out, const glm::mat4& m) {
	out << std::endl;
	out << "(" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << ", " << m[0][3] << "," << std::endl;
	out << " " << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << ", " << m[1][3] << "," << std::endl;
	out << " " << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << ", " << m[2][3] << "," << std::endl;
	out << " " << m[3][0] << ", " << m[3][1] << ", " << m[3][2] << ", " << m[3][3] << ")";
	return out;
}
inline std::ostream& operator<< (std::ostream& out, const glm::quat& q) {
	return out << "(w:" << q.w << ", x:" << q.x << ", y:" << q.y << ", z:" << q.z << ")";
}
// Model::...
inline std::ostream& operator<< (std::ostream& out, const Data::Body& b) {
	return out << "Body: " << b.name << ", transform:" << b.globalTransform;
}
inline std::ostream& operator<< (std::ostream& out, const Data::Joint& j) {
	return out << "Joint: " << j.name << ", childName:" << j.childName << ", parentName: " << j.parentName;
}
inline std::ostream& operator<< (std::ostream& out, const Data::Marker& m) {
	return out << "Marker: " << m.name << ", pos: " << m.globalPosition;
}
inline std::ostream& operator<< (std::ostream& out, const Data::Muscle& m) {
	return out << "Muscle: " << m.name << ", activation: " << m.activation << ", linePoints: " << m.lineSet.size();
}

inline void printRotationMatrix(const glm::mat4& M) {
	glm::vec3 euler = glm::degrees(glm::eulerAngles(glm::quat_cast(M)));
	glm::vec3 pos = M[3];
	PRINT("Location: " << pos << "EulerXYZ: " << euler);
}

//
// find nearest triangle to ray
//
inline bool FindRayHitpointOnMesh(const Ray& ray, const std::vector<glm::vec3>& vertices, std::vector<uint>& indices, const glm::mat4& M, glm::vec3& hitpoint, glm::vec2& baryzentricCoordinates) {
	bool hit = false;
	float tHit = std::numeric_limits<float>::infinity();

	for (uint i = 0; i < indices.size(); i += 3) {
		const uint i0 = indices.at(i);
		const uint i1 = indices.at(i + 1);
		const uint i2 = indices.at(i + 2);

		const glm::vec3 p0 = glm::vec3(glm::vec4(vertices.at(i0), 1) * M);
		const glm::vec3 p1 = glm::vec3(glm::vec4(vertices.at(i1), 1) * M);
		const glm::vec3 p2 = glm::vec3(glm::vec4(vertices.at(i2), 1) * M);

		glm::vec2 baryPosition;
		float distance;
		if (glm::intersectRayTriangle(ray.origin, ray.direction, p0, p1, p2, baryPosition, distance)) {
			// get nearest intersection point
			if (distance < tHit) {
				tHit = distance;
				hitpoint = ray.origin + ray.direction * tHit;
				baryzentricCoordinates = baryPosition;
				hit = true;
			}
		}
	}

	return hit;
}


//
// compute vertex-normals for given mesh (vertices+indices)
//
inline std::vector<glm::vec3> computeVertexNormals(const std::vector<glm::vec3>& vertices, const std::vector<uint>& indices) {
	if (indices.size() % 3 != 0) {
		PRINT_ERR("Only triangles are supported.");
	}
	std::vector<glm::vec3> normals(vertices.size(), glm::vec3(0.0));
	for (uint i = 0; i < indices.size(); i += 3) {
		uint i0 = indices.at(i);
		uint i1 = indices.at(i + 1);
		uint i2 = indices.at(i + 2);

		const glm::vec3& v0 = vertices.at(i0);
		const glm::vec3& v1 = vertices.at(i1);
		const glm::vec3& v2 = vertices.at(i2);

		glm::vec3 e0 = v1 - v0;
		glm::vec3 e1 = v2 - v0;
		glm::vec3 faceNormal = glm::cross(e0, e1);

		normals.at(i0) += faceNormal;
		normals.at(i1) += faceNormal;
		normals.at(i2) += faceNormal;
	}

	// normalize the sum to receive the vertex normals
	for (uint i = 0; i < normals.size(); i++) {
		normals.at(i) = glm::normalize(normals.at(i));
	}
	return normals;
}