#include "modelLoader.h"


//
// ModelLoader
//
ModelLoader::ModelLoader(const std::string& filepath) {
	// load collada file containing the model data
	LoadModelData(filepath);
}

ModelLoader::~ModelLoader() {
	SAFE_DELETE(rootNode); // the destructor of the rootNode handles the memory deleltion of its children
}

// initializes the node hirarchy from the assimp import
void ModelLoader::InitNodeHeirarchy(BoneNode* parentNode, const aiNode* n) {
	// for current node: create child nodes for all children
	for (int i = 0; i < n->mNumChildren; i++) {
		const aiNode* childNode = n->mChildren[i];

		std::string boneName = childNode->mName.C_Str();
		glm::mat4 localTransformation = toGLM(childNode->mTransformation);

		BoneNode* node = new BoneNode(
			boneName,
			localTransformation,
			childNode->mNumChildren,
			parentNode
		);
		parentNode->children.push_back(node);

		boneNodeMap[boneName] = node;
		InitNodeHeirarchy(node, n->mChildren[i]);
	}
}


// computes the world-space transformation per bone
// by traversing the boneNode hierarchy to compute the transformations for each bone given a Animation matrix for each bone
std::vector<glm::mat4> ModelLoader::GetBoneTransformations(const std::map<std::string, glm::mat4>& transformationsInBoneSpace) {
	std::vector<glm::mat4> transformPerBone_World(numBones, glm::identity<glm::mat4>());
	GetTransformations_Heirarchy(transformationsInBoneSpace, rootNode, glm::identity<glm::mat4>(), transformPerBone_World);
	return transformPerBone_World;
}
// traverse the boneNode hierarchy to apply the transformations to each bone
// this is how it is done usually for model animations
void ModelLoader::GetTransformations_Heirarchy(const std::map<std::string, glm::mat4>& transformationsInBoneSpace, const BoneNode* node, const glm::mat4& parentTransform, std::vector<glm::mat4>& transformPerBone_World) {
	// B = Root-1 * Parent(..) * Node * Anim * Offset
	glm::mat4 nodeTransformation = node->localTransformation;

	// Animation: transform nodeTransformation here in bone space
	if (FindKeyInMap(node->name, transformationsInBoneSpace)) {
		nodeTransformation *= transformationsInBoneSpace.at(node->name);
	}
	glm::mat4 globalTransformation = parentTransform * nodeTransformation;

	if (FindKeyInMap(node->name, boneMapping)) {
		int boneID = boneMapping.at(node->name);
		glm::mat4 MeshToBone = offsetsPerBone.at(boneID);
		transformPerBone_World.at(boneID) = glm::inverse(rootNode->localTransformation) * globalTransformation * MeshToBone;
	}
	for (int i = 0; i < node->numChildren; i++) {
		GetTransformations_Heirarchy(transformationsInBoneSpace, node->children[i], globalTransformation, transformPerBone_World);
	}
}


// call with xyz vectors for scale
std::vector<glm::mat4> ModelLoader::GetBoneTransformations_scalePerBone(const std::map<std::string, glm::vec3>& scalesPerBone) {
	std::map<std::string, glm::mat4> scales;
	for (const auto& s : scalesPerBone) {
		scales[s.first] = glm::scale(s.second);
	}
	return GetBoneTransformations_scalePerBone(scales);
}
// call with scaling matrix
std::vector<glm::mat4> ModelLoader::GetBoneTransformations_scalePerBone(const std::map<std::string, glm::mat4>& scalesPerBone) {
	std::vector<glm::mat4> transformPerBone_World(numBones, glm::identity<glm::mat4>());
	GetTransformations_Heirarchy_scalePerBone(scalesPerBone, rootNode, glm::identity<glm::mat4>(), glm::identity<glm::mat4>(), transformPerBone_World);
	return transformPerBone_World;
}
// apply scale to each bone individually by traversing the hierarchy
// note that the bone positions of child-bones are re-calculated when a parent bone is scaled
void ModelLoader::GetTransformations_Heirarchy_scalePerBone(const std::map<std::string, glm::mat4>& scalesPerBone, const BoneNode* node, const glm::mat4& parentTransform, const glm::mat4& parentScale, std::vector<glm::mat4>& transformPerBone_World) {
	// B = Root-1 * Parent(..) * Node * BoneScale * Offset * ParentScale-1
	glm::mat4 nodeTransformation = node->localTransformation;
	glm::mat4 boneScale = glm::identity<glm::mat4>();

	// Animation: transform nodeTransformation here in bone space
	if (FindKeyInMap(node->name, scalesPerBone)) {
		boneScale = scalesPerBone.at(node->name);
		nodeTransformation *= boneScale;
	}
	glm::mat4 globalTransformation = parentTransform * nodeTransformation * glm::inverse(parentScale);

	if (FindKeyInMap(node->name, boneMapping)) {
		int boneID = boneMapping.at(node->name);
		glm::mat4 MeshToBone = offsetsPerBone.at(boneID);
		transformPerBone_World.at(boneID) = glm::inverse(rootNode->localTransformation) * globalTransformation * MeshToBone;
	}
	for (int i = 0; i < node->numChildren; i++) {
		GetTransformations_Heirarchy_scalePerBone(scalesPerBone, node->children[i], globalTransformation, boneScale, transformPerBone_World);
	}
}


// loads all neccessary data from the model file
void ModelLoader::LoadModelData(const std::string& filepath) {
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(filepath, aiProcess_FixInfacingNormals);
	if (!scene) {
		PRINT_ERR(importer.GetErrorString());
	}
	if (scene->mNumMeshes != 1) {
		PRINT_ERR("Only one mesh is supported.");
	}

	aiMesh* m = scene->mMeshes[0];
	if (!m->HasPositions() || !m->HasFaces()) {
		PRINT_ERR("The mesh does not have any vertices or faces.");
	}

	// vertices
	vertices.clear();
	normals.clear();
	for (uint i = 0; i < m->mNumVertices; ++i) {
		const aiVector3D v = m->mVertices[i];
		vertices.push_back(toGLM(v));

		// normals
		if (m->HasNormals()) {
			const aiVector3D n = m->mNormals[i];
			normals.push_back(glm::normalize(toGLM(n)));
		}
	}
	// indices
	indices.clear();
	for (uint i = 0; i < m->mNumFaces; ++i) {
		aiFace f = m->mFaces[i];
		if (f.mNumIndices != 3) {
			PRINT_ERR("Only triangles are supported.");
		}
		indices.push_back(f.mIndices[0]);
		indices.push_back(f.mIndices[1]);
		indices.push_back(f.mIndices[2]);
	}


	// bones
	numBones = m->mNumBones;
	boneWeights.resize(vertices.size());
	offsetsPerBone.resize(numBones);
	boneList.resize(numBones);
	for (uint boneID = 0; boneID < m->mNumBones; boneID++) {
		const aiBone* b = m->mBones[boneID];
		const std::string& boneName = b->mName.C_Str();
		const aiMatrix4x4& M = b->mOffsetMatrix;

		offsetsPerBone.at(boneID) = toGLM(M);
		boneMapping[boneName] = boneID;
		boneList.at(boneID) = boneName;

		if (b->mNumWeights > 0) {
			for (uint j = 0; j < b->mNumWeights; j++) {
				uint vertexID = b->mWeights[j].mVertexId;
				float weight = b->mWeights[j].mWeight;

				boneWeights.at(vertexID).AddBoneWeight(boneID, weight);
			}
		} else {
			PRINT_WAR("Bone " + boneName + " doesn't have any weights.");
		}
	}

	in_boneIDs.clear();
	in_weights.clear();
	for (uint i = 0; i < boneWeights.size(); i++) {
		in_boneIDs.push_back(glm::ivec4(0));
		in_weights.push_back(glm::vec4(0.0));

		for (uint j = 0; j < 4; j++) {
			in_boneIDs.at(i)[j] = boneWeights.at(i).boneIDs[j];
			in_weights.at(i)[j] = boneWeights.at(i).weights[j];
		}
	}

	// bone nodes
	rootNode = new BoneNode(
		scene->mRootNode->mName.C_Str(),
		toGLM(scene->mRootNode->mTransformation),
		scene->mRootNode->mNumChildren,
		NULL
	);
	InitNodeHeirarchy(rootNode, scene->mRootNode);
}

