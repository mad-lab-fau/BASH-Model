#pragma once

#include "util.h"

#include <OpenSim/OpenSim.h>

//
// conversions: SimTK -> GLM
//
inline glm::vec3 toGLM(const SimTK::Vec3& I) {
	return glm::vec3(I[0], I[1], I[2]);
}
inline glm::mat3 toGLM(const SimTK::Mat33& I) {
	glm::mat3 O = glm::identity<glm::mat3>();
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			O[x][y] = I[y][x]; // swap row & col
		}
	}
	return O;
}
inline glm::mat4 toGLM(const SimTK::Mat44& I) {
	glm::mat4 O = glm::identity<glm::mat4>();
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			O[x][y] = I[y][x]; // swap row & col
		}
	}
	return O;
}
inline glm::mat4 toGLM(const SimTK::Transform& I) {
	return toGLM(I.toMat44());
}
inline glm::mat4 toGLM(const SimTK::Rotation& I) {
	return toGLM(I.asMat33());
}

//
// conversions: GLM -> SimTK
//
inline SimTK::Vec3 toSimTK(const glm::vec3& I) {
	return SimTK::Vec3(I[0], I[1], I[2]);
}
inline SimTK::Mat33 toSimTK(const glm::mat3& I) {
	SimTK::Mat33 O;
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			O[x][y] = I[y][x]; // swap row & col
		}
	}
	return O;
}
inline SimTK::Mat44 toSimTK(const glm::mat4& I) {
	SimTK::Mat44 O;
	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			O[x][y] = I[y][x]; // swap row & col
		}
	}
	return O;
}



//
// OSIM
//
class OSIM {
private:
	bool init = false;
	OpenSim::Model osimModel;
	SimTK::State state;
	OpenSim::Storage motStorage;
	std::vector<std::map<std::string, double>> muscleActivations;
	std::map<std::string, glm::vec3> scaleFactors;

	OpenSim::Storage LoadStorageFile(const std::string& filepath);
	void LoadMuscleActivationFile(const std::string& filepath);
	void SetStateToDefault();
	void SetStateToFrame(int frameID);
	std::map<std::string, Data::Marker> LoadMarkerData();
	std::map<std::string, Data::Body> LoadBodyData();
	std::map<std::string, Data::Joint> LoadJointData();
	std::map<std::string, Data::Muscle> LoadMuscleData();

public:
	OSIM(const std::string& filepath_osim, const std::string& filepath_scale = "", const std::string& filepath_mot = "", const std::string& filepath_sto = "");

	int GetNumFrames();

	void LoadScaleFactors(const std::string& filepath_scaleConfig);
	void InverseKinematics(const std::string& filepath_input_trc, const std::string& filepath_output_mot);

	std::map<std::string, glm::vec3> GetScaleFactors();
	std::map<std::string, Data::Marker> GetMarkers();
	std::map<std::string, Data::Marker> GetMarkers(int frameID);
	std::map<std::string, Data::Body> GetBodies();
	std::map<std::string, Data::Body> GetBodies(int frameID);
	std::map<std::string, Data::Joint> GetJoints();
	std::map<std::string, Data::Joint> GetJoints(int frameID);
	std::map<std::string, Data::Muscle> GetMuscles();
	std::map<std::string, Data::Muscle> GetMuscles(int frameID);
};

