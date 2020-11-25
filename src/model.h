#pragma once

#include "util.h"
#include "settings.h"
#include "modelLoader.h"
#include <render\shader.h>


//
// Model
//
class Model {
private:
	bool init = false;
	Model() {};
	Model(const Model&) = delete;
	Model& operator=(const Model&) = delete;

	void CheckInit();

	void LoadScapeMarkers(const std::string& filepath);
	void ComputeModelPose_Default();
	void ComputeModelPose_Scaled();
	void ComputeModelPose_Transformed();
	void ComputeModelPose_ScapeSpace();

	std::vector<glm::vec3> TransformMeshInPose(const std::vector<glm::vec3>& vertices, const std::vector<glm::mat4>& B);
	std::map<std::string, Data::Marker> TransformMarkersInPose(const std::map<std::string, Data::Marker>& markers, const std::vector<glm::mat4>& B);
	std::vector<glm::mat4> GetBoneTransformations_mapping(int frameID);
	std::map<std::string, glm::vec3> EstimateBoneScales(const std::map<std::string, Data::Marker>& markersOsim, const std::map<std::string, Data::Marker>& markersScape);

	void LoadFixedModelData(const std::string& filepath);
	bool LoadFromCache(const std::string& filepath, std::vector<glm::vec3>& vertices);
	void WriteToCache(const std::string& filepath, const std::vector<glm::vec3>& vertices);

	std::vector<uint> kNN(const std::vector<glm::vec3>& points, const std::vector<glm::vec3>& lineSet);

public:
	static Model& GetInstance() {
		static Model instance;
		return instance;
	}

	ModelLoader* modelData = NULL;
	void InitModel();

	// stats
	int numFrames = 0;
	int numVertices = 0;
	int numIndices = 0;
	int numMuscles = 0;

	// mesh data
	std::vector<glm::vec3> vertices_default;
	std::vector<glm::vec3> vertices_scaled;
	std::vector<std::vector<glm::vec3>> vertices_transformedPerFrame;
	std::vector<std::vector<glm::vec3>> vertices_scapeSpacePerFrame;

	std::vector<std::vector<float>> muscleActivityPerFrame;

	std::map<std::string, Data::Marker> markers_default;
	std::map<std::string, Data::Marker> markers_scaled;
	std::vector<std::map<std::string, Data::Marker>> markers_transformedPerFrame;

	// OSIM data
	std::map<std::string, glm::vec3> scaleFactors;

	std::map<std::string, Data::Body> osimBodies_default;
	std::map<std::string, Data::Joint> osimJoints_default;
	std::map<std::string, Data::Marker> osimMarkers_default;
	std::map<std::string, Data::Muscle> osimMuscles_default;

	std::map<std::string, Data::Body> osimBodies_scapePose;
	std::map<std::string, Data::Joint> osimJoints_scapePose;
	std::map<std::string, Data::Marker> osimMarkers_scapePose;
	std::map<std::string, Data::Muscle> osimMuscles_scapePose;

	std::vector<std::map<std::string, Data::Body>> osimBodies_perFrame;
	std::vector<std::map<std::string, Data::Joint>> osimJoints_perFrame;
	std::vector<std::map<std::string, Data::Marker>> osimMarkers_perFrame;
	std::vector<std::map<std::string, Data::Muscle>> osimMuscles_perFrame;
};

