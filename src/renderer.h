#pragma once

#include "util.h"
#include "settings.h"
#include "io.h"
#include "model.h"

#include <render/mesh.h>

//
// States
//
enum class ModelState {
	Default,
	Scaled,
	Transformed,
	ScapeSpace
};
enum class OsimState {
	Default,
	Scape,
	Animation
};

//
// ModelRenderer
//
class ModelRenderer {
private:
	// mesh container to switch between processing states in real time
	Mesh* mesh_default;
	Mesh* mesh_scaled;
	std::vector<Mesh*> mesh_transformedPerFrame;
	std::vector<Mesh*> mesh_scapeSpacePerFrame;
	std::vector<Mesh*> mesh_muscleColorsPerFrame;

	Mesh* meshMarkers_default;
	Mesh* meshMarkers_scaled;
	std::vector<Mesh*> meshMarkers_transformedPerFrame;

	ModelState modelState = ModelState::ScapeSpace;

public:
	ModelRenderer();
	~ModelRenderer();

	Bounds<glm::vec3> bounds;
	glm::mat4 M = glm::identity<glm::mat4>();

	void DrawMarkers(int frameID);
	void Draw(int frameID);

	void SetModelState(ModelState modelState);
	inline ModelState GetModelState() { return modelState; };
};


//
// OsimRenderer
//
class OsimRenderer {
private:
	Mesh* mesh_bonesDefault;
	Mesh* mesh_jointsDefault;
	Mesh* mesh_markersDefault;
	std::vector<Mesh*> mesh_musclesDefault;

	Mesh* mesh_bonesScapePose;
	Mesh* mesh_jointsScapePose;
	Mesh* mesh_markersScapePose;
	std::vector<Mesh*> mesh_musclesScapePose;

	std::vector<Mesh*> mesh_bonesPerFrame;
	std::vector<Mesh*> mesh_jointsPerFrame;
	std::vector<Mesh*> mesh_markersPerFrame;
	std::vector<std::vector<Mesh*>> mesh_musclesPerFrame;
	std::vector<std::vector<glm::vec3>> muscleActivationColors;

	OsimState osimState = OsimState::Animation;

	void Draw_Default();
	void Draw_Scape();
	void Draw_Animation(int frameID);

public:
	OsimRenderer();
	~OsimRenderer();

	glm::mat4 M = glm::identity<glm::mat4>();

	void Draw(int frameID);

	inline void SetOsimState(OsimState s) { osimState = s; }
	inline OsimState GetOsimState() { return osimState; };
};

