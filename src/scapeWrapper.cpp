#include "scapeWrapper.h"

#include <SCAPE/SimpleMesh.h>
#include <SCAPE/rigid_transform.h>
#include <SCAPE/operators/alignPointClouds.h>
#include <SCAPE/FileIO/MeshWriter.h>

//
// SCAPE
//
void SCAPE::InitSCAPEModel(const std::string& binFile, const std::string& dataDirectory) {
	init = true;

	scapeModel = SCAPE_Model();
	scapeModel.readFromFile(binFile.c_str(), dataDirectory.c_str(), true);
	PRINT("Successfully loaded SCAPE model from: " << dataDirectory);
}

void SCAPE::CheckInit() {
	if (!init) {
		PRINT_ERR("SCAPE must be initialized first!");
	}
}


// create a mesh from given SCAPE parameters
SimpleMesh SCAPE::CreateMeshFromSCAPEParams(const std::vector<double>& scapeParameters) {
	CheckInit();
	return SimpleMesh(scapeModel.reconstructFromSCAPEParams(scapeParameters), scapeModel.topologically_correct_mesh_triangles);
}

// project to SCAPE space
std::vector<vec3d> SCAPE::ProjectToScapeSpace(const std::vector<glm::vec3>& vertices) {
	CheckInit();

	// custom shape params
	std::vector<double> shapeParams(NumShapeParams(), 0.0);

	//return scapeModel.projectToSCAPESpace(toVec3d(vertices), false, shapeParams);
	return scapeModel.projectToSCAPESpace(toVec3d(vertices), true); // pose only
}


// project a input mesh to the SCAPE space
// NOTE: indices stay always the same!
std::vector<glm::vec3> SCAPE::FitScapeModel(const std::vector<glm::vec3>& vertices) {
	CheckInit();

	// input mesh
	SimpleMesh inputMesh = SimpleMesh(vertices, scapeModel.topologically_correct_mesh_triangles);

	// project to scape space (pose only!)
	//SimpleMesh transformedMesh = SimpleMesh(ProjectToScapeSpace_Pose(vertices), scapeModel.topologically_correct_mesh_triangles);
	SimpleMesh transformedMesh = SimpleMesh(ProjectToScapeSpace(vertices), scapeModel.topologically_correct_mesh_triangles);

	// fix rigid alignment (Scale+Rotation+Translation) because SCAPE has only relative joint rotations and it's right foot is fixed in the origin
	RigidTransformWithScale<vec3d> trans = getOptimalAlignment(&inputMesh, &transformedMesh);
	transformMesh(&transformedMesh, trans);

	return toGLM(transformedMesh);
}


