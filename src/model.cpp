#include "model.h"

#include "io.h"
#include "scapeWrapper.h"
#include "osim.h"

#include <glm/gtx/closest_point.hpp>
#include <filesystem>


//
// Model
//
void Model::InitModel() {
	init = true;

	// load model mesh data and fix assimp vertex order
	modelData = new ModelLoader(Settings::GetInstance().baselineModelDir + FILENAME_BASELINEMODEL);
	LoadFixedModelData(Settings::GetInstance().baselineModelDir);

	// load input model
	OSIM osimInput(Settings::GetInstance().inputModelOSIM, Settings::GetInstance().inputModelScale, Settings::GetInstance().inputModelMOT, Settings::GetInstance().inputModelSTO);
	scaleFactors = osimInput.GetScaleFactors();

	// load the SCAPE model
	SCAPE::GetInstance().InitSCAPEModel(std::string(FILEPATH_SCAPE_DATADIRECTORY) + std::string(FILENAME_SCAPE_BINARYDATA), FILEPATH_SCAPE_DATADIRECTORY);

	// OSIM model in default pose
	OSIM osimZeros(Settings::GetInstance().inputModelOSIM);
	osimBodies_default = osimZeros.GetBodies();
	osimJoints_default = osimZeros.GetJoints();
	osimMarkers_default = osimZeros.GetMarkers();
	osimMuscles_default = osimZeros.GetMuscles();

	// load SCAPE mesh and markers in the default pose
	ComputeModelPose_Default();

	// scale model (mesh and markers)
	ComputeModelPose_Scaled();

	// write SCAPE markers to a .trc file that will be input for IK
	IO::GetInstance().WriteTRC_Markers(std::string(FILEPATH_CACHE_MAPPING) + std::string(FILENAME_MARKERSONSCAPE_TRC), markers_scaled);

	// compute inverse kinematics and store the inverse pose mapping matrix in a .mot file
	osimInput.InverseKinematics(std::string(FILEPATH_CACHE_MAPPING) + std::string(FILENAME_MARKERSONSCAPE_TRC), std::string(FILEPATH_CACHE_MAPPING) + std::string(FILENAME_POSEMAPPING_MOT));

	// receive OSIM model that stores the inverse pose mapping transformed components in the first frame
	OSIM osimInScapePose(Settings::GetInstance().inputModelOSIM, Settings::GetInstance().inputModelScale, std::string(FILEPATH_CACHE_MAPPING) + std::string(FILENAME_POSEMAPPING_MOT));
	osimBodies_scapePose = osimInScapePose.GetBodies(0);
	osimJoints_scapePose = osimInScapePose.GetJoints(0);
	osimMarkers_scapePose = osimInScapePose.GetMarkers(0);
	osimMuscles_scapePose = osimInScapePose.GetMuscles();

	// store input model data in global storage (for each frame)
	for (int frameID = 0; frameID < osimInput.GetNumFrames(); frameID++) {
		// OSIM
		std::map<std::string, Data::Body> osim_bodies = osimInput.GetBodies(frameID);
		std::map<std::string, Data::Joint> osim_joints = osimInput.GetJoints(frameID);
		std::map<std::string, Data::Marker> osim_markers = osimInput.GetMarkers(frameID);
		std::map<std::string, Data::Muscle> osim_muscles = osimInput.GetMuscles(frameID);

		numFrames++;
		osimBodies_perFrame.push_back(osim_bodies);
		osimJoints_perFrame.push_back(osim_joints);
		osimMarkers_perFrame.push_back(osim_markers);
		osimMuscles_perFrame.push_back(osim_muscles);

		// DEBUG
		if (frameID >= 20) {
			//break;
		}
	}
	// muscle IDs
	numMuscles = osimMuscles_perFrame.at(0).size();

	// compute transformed vertices
	ComputeModelPose_Transformed();

	// project mesh to scape space
	ComputeModelPose_ScapeSpace();

	// find area of influence of muscles in scape default pose
	PRINT("Muscle: Find area of influence on surface");
	std::map<std::string, std::vector<uint>> muscleInfluencingVertices;
	for (const auto& m : osimMuscles_scapePose) {
		const Data::Muscle& muscle = m.second;

		if (MUSCLE_MODE == MuscleMode::LoadFromFile) {
			// load from file	
			std::string filepath = std::string(FILEPATH_MUSCLES) + muscle.name + ".txt";
			if (FileExists(filepath)) {
				muscleInfluencingVertices[muscle.name] = IO::GetInstance().ReadVectorDataFromFile<uint>(filepath);
			} else {
				PRINT_WAR("File does not exist: " + filepath);
			}
		} else if (MUSCLE_MODE == MuscleMode::kNN) {
			// kNN for every muscle
			muscleInfluencingVertices[muscle.name] = kNN(vertices_scaled, muscle.lineSet);
		}
	}

	// compute per-vertex colors to visualize muscle activity on the mesh
	PRINT("Muscle: Compute color values to visualize musclular acitvity");
	muscleActivityPerFrame.resize(numFrames);
	float maxMuscleActivity = 0.0;
	for (int frameID = 0; frameID < numFrames; frameID++) {
		std::vector<float> muscleActivities(numVertices, 0.0);

		if (Settings::GetInstance().visualizeMuscleActivity) {
			// for every muscle...
			for (const auto& m : osimMuscles_perFrame.at(frameID)) {
				const Data::Muscle& muscle = m.second;
				if (ONLY_RENDER_SINGLE_MUSCLE != "" && muscle.name != ONLY_RENDER_SINGLE_MUSCLE) continue; // debugging

				if (FindKeyInMap(muscle.name, muscleInfluencingVertices)) {
					for (uint vertexID : muscleInfluencingVertices.at(muscle.name)) {
						muscleActivities.at(vertexID) += muscle.activation * muscle.maxForce;

						maxMuscleActivity = glm::max(maxMuscleActivity, muscleActivities.at(vertexID));
					}
				}
			}

			// fit to color range [0..1] in regard of maximal possible force
			for (auto& muscleActivity : muscleActivities) {
				muscleActivity /= MUSCLE_MAX_ACTIVITY;
			}
		}

		muscleActivityPerFrame.at(frameID) = muscleActivities;
	}
}

std::vector<uint> Model::kNN(const std::vector<glm::vec3>& points, const std::vector<glm::vec3>& lineSet) {
	std::vector<uint> pointsHit;

	float minDistance = 1000.0;
	float maxDistance = 0.0;

	// iterate through all points
	for (uint pointID = 0; pointID < points.size(); pointID++) {
		const glm::vec3& vertex = points.at(pointID);

		// check for every line segment of the muscle
		for (uint i = 0; i < lineSet.size() - 1; i++) {
			Line line;
			line.A = lineSet.at(i);
			line.B = lineSet.at(i + 1);

			glm::vec3 pointOnLine = glm::closestPointOnLine(vertex, line.A, line.B);
			float distance = glm::distance(vertex, pointOnLine);

			if (distance < MUSCLE_MAX_DISTANCE_TO_SURFACE) { // TODO: encorporate scale...?
				pointsHit.push_back(pointID);
			}

			minDistance = glm::min(minDistance, distance);
			maxDistance = glm::max(maxDistance, distance);
		}
	}
	return pointsHit;
}

void Model::CheckInit() {
	if (!init) {
		PRINT_ERR("Model must be initialized first!");
	}
}

void Model::LoadScapeMarkers(const std::string& filepath) {
	CheckInit();

	markers_default = IO::GetInstance().ReadOBJ_Markers(filepath);

	// attach SCAPE markers to model bones
	for (auto& m : markers_default) {
		Data::Marker& marker = m.second;

		if (!FindKeyInMap(marker.name, osimMarkers_default)) {
			PRINT_ERR("No corresponding marker was found: " + marker.name);
		}
		// get informations from corresponding OSIM marker and compute local position respective to their bone
		const Data::Marker& markerOsim = osimMarkers_default.at(marker.name);
		const std::string& boneName = markerOsim.attachedBody;
		int boneID = modelData->boneMapping.at(boneName);
		glm::mat4 MeshToBone = modelData->offsetsPerBone.at(boneID);

		// complete missing attributes
		marker.attachedBody = boneName;
		marker.localPosition = glm::vec3(glm::vec4(marker.globalPosition, 1.0) * MeshToBone);
	}
}

void Model::ComputeModelPose_Default() {
	CheckInit();
	PRINT("Retrieving model's default pose");

	vertices_default = modelData->vertices;
	LoadScapeMarkers(Settings::GetInstance().baselineModelDir + FILENAME_BASELINEMODEL_MARKERS);

	if (CACHE_ALL_MESHES) {
		// just write the mesh
		std::string fileCachedMesh = Settings::GetInstance().filepathModelCache + std::string(MESH_NAME_DEFAULT) + ".ply";
		WriteToCache(fileCachedMesh, vertices_default);
	}
}

void Model::ComputeModelPose_Scaled() {
	CheckInit();
	PRINT("Scaling model");

	for (const auto& s : scaleFactors) {
		PRINT("Scale factor " << s.first << ": " << s.second);
	}

	// optionally, estimate the scalings per bone from marker data and override existing
	if (ESTIMATE_SCALEFACTORS_MARKERDATA) {
		PRINT("Estimating scale factors based on marker data. Overriding existing values!");
		scaleFactors = EstimateBoneScales(osimMarkers_default, markers_default);
	}

	std::vector<glm::mat4> S = modelData->GetBoneTransformations_scalePerBone(scaleFactors);
	vertices_scaled = TransformMeshInPose(vertices_default, S);
	markers_scaled = TransformMarkersInPose(markers_default, S);

	if (CACHE_ALL_MESHES) {
		// just write the mesh
		std::string fileCachedMesh = Settings::GetInstance().filepathModelCache + std::string(MESH_NAME_SCALED) + ".ply";
		WriteToCache(fileCachedMesh, vertices_scaled);
	}
}

void Model::ComputeModelPose_Transformed() {
	CheckInit();

	vertices_transformedPerFrame.resize(numFrames);
	markers_transformedPerFrame.resize(numFrames);
	for (int frameID = 0; frameID < numFrames; frameID++) {
		PRINT("Transforming model's pose for frame " << frameID);
		std::vector<glm::mat4> B = GetBoneTransformations_mapping(frameID);
		markers_transformedPerFrame.at(frameID) = TransformMarkersInPose(markers_scaled, B);

		if (CACHE_ALL_MESHES) {
			// load or compute & store from cache
			std::string fileCachedMesh = Settings::GetInstance().filepathModelCache + std::string(MESH_NAME_TRANSFORMED) + "_" + to_zero_lead(frameID, 4) + ".ply";
			if (!LoadFromCache(fileCachedMesh, vertices_transformedPerFrame.at(frameID))) {
				vertices_transformedPerFrame.at(frameID) = TransformMeshInPose(vertices_scaled, B);
				WriteToCache(fileCachedMesh, vertices_transformedPerFrame.at(frameID));
			}
		} else {
			vertices_transformedPerFrame.at(frameID) = TransformMeshInPose(vertices_scaled, B);
		}
	}
}

void Model::ComputeModelPose_ScapeSpace() {
	CheckInit();

	vertices_scapeSpacePerFrame.resize(numFrames);
	for (int frameID = 0; frameID < numFrames; frameID++) {
		PRINT("Projecting model to ScapeSpace for frame " << frameID);

		if (CACHE_ALL_MESHES) {
			// load or compute & store from cache
			std::string fileCachedMesh = Settings::GetInstance().filepathModelCache + std::string(MESH_NAME_SCAPE) + "_" + to_zero_lead(frameID, 4) + ".ply";
			if (!LoadFromCache(fileCachedMesh, vertices_scapeSpacePerFrame.at(frameID))) {
				vertices_scapeSpacePerFrame.at(frameID) = SCAPE::GetInstance().FitScapeModel(vertices_transformedPerFrame.at(frameID));
				WriteToCache(fileCachedMesh, vertices_scapeSpacePerFrame.at(frameID));
			}
		} else {
			vertices_scapeSpacePerFrame.at(frameID) = SCAPE::GetInstance().FitScapeModel(vertices_transformedPerFrame.at(frameID));
		}
	}
}

// apply the given bone transformations to the given vertices (this is usually done in the vertex-shader, but we need the transformed vertices here to apply the SCAPE projection)
std::vector<glm::vec3> Model::TransformMeshInPose(const std::vector<glm::vec3>& vertices, const std::vector<glm::mat4>& B) {
	CheckInit();

	std::vector<glm::vec3> transformedVertices(vertices.size());
	for (int vertexID = 0; vertexID < vertices.size(); vertexID++) {
		glm::vec3 in_position = vertices.at(vertexID);
		glm::ivec4 in_boneID = modelData->in_boneIDs.at(vertexID);
		glm::vec4 in_weight = modelData->in_weights.at(vertexID);

		glm::mat4 M_bone = B[in_boneID[0]] * in_weight[0] + B[in_boneID[1]] * in_weight[1] + B[in_boneID[2]] * in_weight[2] + B[in_boneID[3]] * in_weight[3];
		glm::vec4 T_position = M_bone * glm::vec4(in_position, 1.0);

		transformedVertices.at(vertexID) = glm::vec3(T_position);
	}
	return transformedVertices;
}

// apply the give bone transformations to a set of markers
std::map<std::string, Data::Marker> Model::TransformMarkersInPose(const std::map<std::string, Data::Marker>& markers, const std::vector<glm::mat4>& B) {
	CheckInit();

	std::map<std::string, Data::Marker> transformedMarkers;
	for (const auto& m : markers) {
		const Data::Marker& marker = m.second;
		int boneID = modelData->boneMapping.at(marker.attachedBody);

		glm::vec4 T_position = B.at(boneID) * glm::vec4(marker.globalPosition, 1.0);

		Data::Marker transformedMarker;
		transformedMarker.name = marker.name;
		transformedMarker.attachedBody = marker.attachedBody;
		transformedMarker.localPosition = marker.localPosition;
		transformedMarker.globalPosition = glm::vec3(T_position);
		transformedMarkers[marker.name] = transformedMarker;
	}
	return transformedMarkers;
}


// get the transformation matrices for each bone given a frameID
// this method apply the inverse pose mapping matrix beforehand, so SCAPE is in the neutral OSIM pose
// then the OSIM transformation can be applied directly
std::vector<glm::mat4> Model::GetBoneTransformations_mapping(int frameID) {
	CheckInit();

	std::vector<glm::mat4> TransformationPerBone(modelData->numBones, glm::identity<glm::mat4>());

	for (auto& b : osimBodies_perFrame.at(frameID)) {
		const Data::Body& body = b.second;

		if (FindKeyInMap(body.name, modelData->boneMapping)) {
			int boneID = modelData->boneMapping.at(body.name);

			glm::mat4 osimTransformation = body.globalTransform;
			glm::mat4 inversePoseMapping = glm::inverse(osimBodies_scapePose.at(body.name).globalTransform);

			// apply the inverse Pose mapping matrix beforehand, then apply the raw osim transformations
			TransformationPerBone.at(boneID) = osimTransformation * inversePoseMapping;

			//TransformationPerBone.at(boneID) = glm::identity<glm::mat4>(); // default SCAPE pose
		}
	}
	return TransformationPerBone;
}


// estimate scales by comparing OSIM markers to SCAPE markers
std::map<std::string, glm::vec3> Model::EstimateBoneScales(const std::map<std::string, Data::Marker>& markersOsim, const std::map<std::string, Data::Marker>& markersScape) {
	CheckInit();
	std::map<std::string, glm::vec3> bodyScales;

	// TODO: maybe generate automatically? -> all combinations of the markers that belong to a bone
	std::map<std::string, std::vector<std::pair<std::string, std::string>>> measurementSet;
	measurementSet["pelvis"] = { std::make_pair("RASI", "LASI"), std::make_pair("SACR", "RASI"), std::make_pair("SACR", "LASI") };
	measurementSet["torso"] = { std::make_pair("SACR", "C7"), std::make_pair("STRN", "CLAV") };

	measurementSet["humerus_r"] = { std::make_pair("RSHO", "RELB"), std::make_pair("RSHO", "RUPA") };
	measurementSet["ulna_r"] = { std::make_pair("RELB", "RFRA"), std::make_pair("RELB", "RWRB") };
	measurementSet["radius_r"] = { std::make_pair("RELB", "RFRA"), std::make_pair("RELB", "RWRB") };
	measurementSet["hand_r"] = { std::make_pair("RWRB", "RWBA") };

	measurementSet["humerus_l"] = { std::make_pair("LSHO", "LELB"), std::make_pair("LSHO", "LUPA") };
	measurementSet["ulna_l"] = { std::make_pair("LELB", "LFRA"), std::make_pair("LELB", "LWRE") };
	measurementSet["radius_l"] = { std::make_pair("LELB", "LFRA"), std::make_pair("LELB", "LWRE") };
	measurementSet["hand_l"] = { std::make_pair("LWRE", "LWRA") };

	measurementSet["femur_r"] = { std::make_pair("TRO", "KneL"), std::make_pair("KneL", "KneM") };
	measurementSet["tibia_r"] = { std::make_pair("Pat", "AnkM"), std::make_pair("KneL", "KneM"), std::make_pair("AnkL", "AnkM") };
	measurementSet["talus_r"] = { std::make_pair("AnkL", "AnkM"), std::make_pair("AnkL", "Hee") };
	measurementSet["calcn_r"] = { std::make_pair("ToeM", "CalM"), std::make_pair("ToeL", "CalL"), std::make_pair("MFM", "MFL") };
	measurementSet["toes_r"] = { std::make_pair("ToeM", "ToeL"), std::make_pair("ToeM", "Toe2") };

	measurementSet["femur_l"] = { std::make_pair("LTRO", "LKneL"), std::make_pair("LKneL", "LKneM") };
	measurementSet["tibia_l"] = { std::make_pair("LPat", "LAnkM"), std::make_pair("LKneL", "LKneM"), std::make_pair("LAnkL", "LAnkM") };
	measurementSet["talus_l"] = { std::make_pair("LAnkL", "LAnkM"), std::make_pair("LAnkL", "LHee") };
	measurementSet["calcn_l"] = { std::make_pair("LToeM", "LCalM"), std::make_pair("LToeL", "LCalL"), std::make_pair("LMFM", "LMFL") };
	measurementSet["toes_l"] = { std::make_pair("LToeM", "LToeL"), std::make_pair("LToeM", "LToe2") };


	for (const auto& measurement : measurementSet) {
		std::string bodyName = measurement.first;
		const auto& measurements = measurement.second;

		for (const auto& markerPair : measurements) {
			const std::string& markerNameA = markerPair.first;
			const std::string& markerNameB = markerPair.second;
			if (!FindKeyInMap(markerNameA, markersOsim) || !FindKeyInMap(markerNameA, markersScape)) {
				PRINT_ERR("Could not find marker: " + markerNameA);
			}
			if (!FindKeyInMap(markerNameB, markersOsim) || !FindKeyInMap(markerNameB, markersScape)) {
				PRINT_ERR("Could not find marker: " + markerNameB);
			}

			const Data::Marker& markerOsimA = markersOsim.at(markerNameA);
			const Data::Marker& markerScapeA = markersScape.at(markerNameA);
			const Data::Marker& markerOsimB = markersOsim.at(markerNameB);
			const Data::Marker& markerScapeB = markersScape.at(markerNameB);

			float lengthOsim = glm::distance(markerOsimA.globalPosition, markerOsimB.globalPosition);
			float lengthScape = glm::distance(markerScapeA.globalPosition, markerScapeB.globalPosition);

			float scale = lengthOsim / lengthScape;

			// compute uniform scaling
			bodyScales[bodyName] += glm::vec3(scale);
		}

		bodyScales[bodyName] /= measurements.size();

		PRINT("Estimated scale " << bodyName << ": " << bodyScales[bodyName].x);
	}

	return bodyScales;
}



// FIX: assimp reorders vertices...so load the corrected data manually
void Model::LoadFixedModelData(const std::string& filepath) {
	CheckInit();

	modelData->indices.clear();
	modelData->vertices.clear();
	modelData->normals.clear();
	modelData->boneWeights.clear();
	modelData->in_boneIDs.clear();
	modelData->in_weights.clear();

	// mesh
	const Data::Mesh& meshData = IO::GetInstance().ReadOBJ(filepath + "fix_mesh.obj");
	modelData->indices = meshData.indices;
	modelData->vertices = meshData.vertices;
	modelData->normals = meshData.normals;
	numVertices = modelData->vertices.size();
	numIndices = modelData->indices.size();

	// weights
	std::vector<uint> numWeightsPerVertex = IO::GetInstance().ReadVectorDataFromFile<uint>(filepath + "fix_weights_numVertex.txt");
	std::vector<uint> jointsPerWeight = IO::GetInstance().ReadVectorDataFromFile<uint>(filepath + "fix_weights_jointIDs.txt");
	std::vector<float> weights = IO::GetInstance().ReadVectorDataFromFile<float>(filepath + "fix_weights_values.txt");

	int counter = 0;
	for (uint i = 0; i < numWeightsPerVertex.size(); i++) {
		uint numWeights = numWeightsPerVertex.at(i);

		modelData->in_boneIDs.push_back(glm::ivec4(0));
		modelData->in_weights.push_back(glm::vec4(0.0));
		for (int j = 0; j < numWeights; j++) {
			uint boneID = jointsPerWeight.at((counter + j) * 2);
			uint weightID = jointsPerWeight.at((counter + j) * 2 + 1);
			float weight = weights.at(weightID);

			modelData->in_boneIDs.at(i)[j] = boneID;
			modelData->in_weights.at(i)[j] = weight;
		}
		counter += numWeights;
	}

	// bodypart colors: caluculate automatically from boneWeights
	modelData->colors.clear();
	for (uint i = 0; i < modelData->vertices.size(); i++) {
		uint dominantPartID = 0;
		float weight = 0;
		for (uint j = 0; j < 4; j++) {
			if (modelData->in_weights.at(i)[j] > weight) {
				dominantPartID = modelData->in_boneIDs.at(i)[j];
			}
		}
		modelData->colors.push_back(Settings::GetInstance().COLOR_MAP.at(dominantPartID));
	}

	// compute vertex normals
	modelData->normals = computeVertexNormals(modelData->vertices, modelData->indices);
}

bool Model::LoadFromCache(const std::string& filepath, std::vector<glm::vec3>& vertices) {
	if (FileExists(filepath)) {
		Data::Mesh meshData = IO::GetInstance().ReadPLY(filepath);
		vertices = meshData.vertices;
		return true;
	} else {
		return false;
	}
}

void Model::WriteToCache(const std::string& filepath, const std::vector<glm::vec3>& vertices) {
	Data::Mesh meshData;
	meshData.vertices = vertices;
	meshData.indices = modelData->indices; // copy indices
	IO::GetInstance().WritePLY(filepath, meshData);
}

