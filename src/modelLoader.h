#pragma once

#include "util.h"
#include "settings.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

//
// conversions: assimp -> GLM
//
inline glm::vec3 toGLM(const aiVector3D& I) {
	return glm::vec3(I.x, I.y, I.z);
}
inline glm::vec3 toGLM(const aiColor3D& I) {
	return glm::vec3(I.r, I.g, I.b);
}
inline glm::mat3 toGLM(const aiMatrix3x3& I) {
	glm::mat3 O;
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			O[x][y] = I[y][x]; // rows and cols are swapped
		}
	}
	return O;
}
inline glm::mat4 toGLM(const aiMatrix4x4& I) {
	glm::mat4 O;
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			O[x][y] = I[y][x]; // rows and cols are swapped
		}
	}
	return O;
}

//
// conversions: GLM -> assimp
//
inline aiVector3D toAssimp(const glm::vec3& I) {
	return aiVector3D(I.x, I.y, I.z);
}
inline aiMatrix3x3 toAssimp(const glm::mat3& I) {
	aiMatrix3x3 O;
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			O[x][y] = I[y][x]; // rows and cols are swapped
		}
	}
	return O;
}
inline aiMatrix4x4 toAssimp(const glm::mat4& I) {
	aiMatrix4x4 O;
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			O[x][y] = I[y][x]; // rows and cols are swapped
		}
	}
	return O;
}


//
// BoneNode
//
struct BoneNode {
	std::string name;
	glm::mat4 localTransformation; // transformation relative to parent boneNode
	uint numChildren;
	std::vector<BoneNode*> children;
	BoneNode* parent;

	BoneNode(const std::string& name, const glm::mat4& localTransformation, uint numChildren, BoneNode* parent)
		: name(name), localTransformation(localTransformation), numChildren(numChildren), parent(parent) {
		children.clear();
	}
	~BoneNode() {
		for (auto& c : children) {
			SAFE_DELETE(c);
		}
	}
};

//
// BoneWeightsPerVertex
//
struct BoneWeightsPerVertex {
	int MAX_BONES_PER_VERTEX = 4;
	std::vector<uint> boneIDs;
	std::vector<float> weights;

	BoneWeightsPerVertex() {
		boneIDs = std::vector<unsigned int>(MAX_BONES_PER_VERTEX, 0);
		weights = std::vector<float>(MAX_BONES_PER_VERTEX, 0.0f);
	};

	inline void AddBoneWeight(uint boneID, float weight) {
		for (unsigned int i = 0; i < MAX_BONES_PER_VERTEX; i++) {
			if (weights[i] == 0.0) {
				boneIDs[i] = boneID;
				weights[i] = weight;
				return;
			}
		}
		PRINT_ERR("Only " + std::to_string(MAX_BONES_PER_VERTEX) + " Bones per Vertex are allowed.");
	}
};

//
// ModelLoader
//
class ModelLoader {
private:
	void LoadModelData(const std::string& filepath);
	void InitNodeHeirarchy(BoneNode* parentNode, const aiNode* n);

	void GetTransformations_Heirarchy(const std::map<std::string, glm::mat4>& transformationsInBoneSpace, const BoneNode* node, const glm::mat4& parentTransform, std::vector<glm::mat4>& B);
	void GetTransformations_Heirarchy_scalePerBone(const std::map<std::string, glm::mat4>& transformationsInBoneSpace, const BoneNode* node, const glm::mat4& parentTransform, const glm::mat4& parentScale, std::vector<glm::mat4>& B);

public:
	ModelLoader(const std::string& filepath);
	~ModelLoader();

	// mesh data
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec3> colors;
	std::vector<uint> indices;

	// bone weights
	int numBones = 0;
	std::vector<BoneWeightsPerVertex> boneWeights;
	std::vector<glm::ivec4> in_boneIDs;
	std::vector<glm::vec4> in_weights;
	std::vector<glm::mat4> offsetsPerBone; // transforms from mesh space to bone space

	// node hierarchy
	BoneNode* rootNode;
	std::vector<std::string> boneList;
	std::map<std::string, int> boneMapping;
	std::map<std::string, BoneNode*> boneNodeMap;

	// functions
	std::vector<glm::mat4> GetBoneTransformations(const std::map<std::string, glm::mat4>& transformationsInBoneSpace);
	std::vector<glm::mat4> GetBoneTransformations_scalePerBone(const std::map<std::string, glm::vec3>& scalesPerBone);
	std::vector<glm::mat4> GetBoneTransformations_scalePerBone(const std::map<std::string, glm::mat4>& scalesPerBone);
};


