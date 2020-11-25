#include "renderer.h"

#include <render/camera.h>
#include <render/shader.h>


//
// ModelRenderer
//
ModelRenderer::ModelRenderer() {
	Model* model = &Model::GetInstance();
	ModelLoader* modelData = model->modelData;

	//
	// markers
	//
	meshMarkers_default = Mesh::CreatePoints(toPointVector(model->markers_default));
	meshMarkers_scaled = Mesh::CreatePoints(toPointVector(model->markers_scaled));
	for (int frameID = 0; frameID < model->vertices_transformedPerFrame.size(); frameID++) {
		Mesh* mesh_t = Mesh::CreatePoints(toPointVector(model->markers_transformedPerFrame.at(frameID)));
		meshMarkers_transformedPerFrame.push_back(mesh_t);
	}

	//
	// mesh - default
	//
	mesh_default = new Mesh("MeshDefault");
	mesh_default->AddBuffer("in_position", static_cast<uint>(model->vertices_default.size()), 3, GL_FLOAT, GL_STATIC_DRAW, model->vertices_default.data());
	mesh_default->AddBuffer("in_normal", static_cast<uint>(modelData->normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->normals.data());
	mesh_default->AddBuffer("in_boneColor", static_cast<uint>(modelData->colors.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->colors.data());
	mesh_default->AddBuffer("in_muscleActivity", static_cast<uint>(model->numVertices), 1, GL_FLOAT, GL_STATIC_DRAW, std::vector<float>(model->numVertices, 0.0).data());
	mesh_default->AddIndexBuffer(static_cast<uint>(modelData->indices.size()), GL_STATIC_DRAW, modelData->indices.data());
	mesh_default->bounds.ExpandBy(model->vertices_default);

	//
	// mesh - scaled
	//
	mesh_scaled = new Mesh("MeshScaled");
	mesh_scaled->AddBuffer("in_position", static_cast<uint>(model->vertices_scaled.size()), 3, GL_FLOAT, GL_STATIC_DRAW, model->vertices_scaled.data());
	mesh_scaled->AddBuffer("in_normal", static_cast<uint>(modelData->normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->normals.data());
	mesh_scaled->AddBuffer("in_boneColor", static_cast<uint>(modelData->colors.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->colors.data());
	mesh_scaled->AddBuffer("in_muscleActivity", static_cast<uint>(model->numVertices), 1, GL_FLOAT, GL_STATIC_DRAW, std::vector<float>(model->numVertices, 0.0).data());
	mesh_scaled->AddIndexBuffer(static_cast<uint>(modelData->indices.size()), GL_STATIC_DRAW, modelData->indices.data());
	mesh_scaled->bounds.ExpandBy(model->vertices_scaled);

	//
	// mesh - transformed
	//
	for (int frameID = 0; frameID < model->vertices_transformedPerFrame.size(); frameID++) {
		Mesh* mesh_t = new Mesh("MeshScapeSpace");
		mesh_t->AddBuffer("in_position", static_cast<uint>(model->vertices_transformedPerFrame.at(frameID).size()), 3, GL_FLOAT, GL_STATIC_DRAW, model->vertices_transformedPerFrame.at(frameID).data());
		mesh_t->AddBuffer("in_normal", static_cast<uint>(modelData->normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->normals.data());
		mesh_t->AddBuffer("in_boneColor", static_cast<uint>(modelData->colors.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->colors.data());
		mesh_t->AddBuffer("in_muscleActivity", static_cast<uint>(model->muscleActivityPerFrame.at(frameID).size()), 1, GL_FLOAT, GL_STATIC_DRAW, model->muscleActivityPerFrame.at(frameID).data());
		mesh_t->AddIndexBuffer(static_cast<uint>(modelData->indices.size()), GL_STATIC_DRAW, modelData->indices.data());
		mesh_t->bounds.ExpandBy(model->vertices_transformedPerFrame.at(frameID));
		mesh_transformedPerFrame.push_back(mesh_t);
	}

	//
	// mesh - scapeSpace
	//
	for (int frameID = 0; frameID < model->vertices_scapeSpacePerFrame.size(); frameID++) {
		Mesh* mesh_s = new Mesh("MeshScapeSpace");
		mesh_s->AddBuffer("in_position", static_cast<uint>(model->vertices_scapeSpacePerFrame.at(frameID).size()), 3, GL_FLOAT, GL_STATIC_DRAW, model->vertices_scapeSpacePerFrame.at(frameID).data());
		mesh_s->AddBuffer("in_normal", static_cast<uint>(modelData->normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->normals.data());
		mesh_s->AddBuffer("in_boneColor", static_cast<uint>(modelData->colors.size()), 3, GL_FLOAT, GL_STATIC_DRAW, modelData->colors.data());
		mesh_s->AddBuffer("in_muscleActivity", static_cast<uint>(model->muscleActivityPerFrame.at(frameID).size()), 1, GL_FLOAT, GL_STATIC_DRAW, model->muscleActivityPerFrame.at(frameID).data());
		mesh_s->AddIndexBuffer(static_cast<uint>(modelData->indices.size()), GL_STATIC_DRAW, modelData->indices.data());
		mesh_s->bounds.ExpandBy(model->vertices_scapeSpacePerFrame.at(frameID));
		mesh_scapeSpacePerFrame.push_back(mesh_s);
	}

	SetModelState(ModelState::ScapeSpace);

}

ModelRenderer::~ModelRenderer() {
	SAFE_DELETE(mesh_default);
	SAFE_DELETE(mesh_scaled);
	for (auto& mesh : mesh_transformedPerFrame) {
		SAFE_DELETE(mesh);
	}
	for (auto& mesh : mesh_scapeSpacePerFrame) {
		SAFE_DELETE(mesh);
	}
}


void ModelRenderer::DrawMarkers(int frameID) {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("color");

	shader->Use();
	{
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);

		// markers
		glPointSize(6.0f);
		shader->Uniform("color", glm::vec3(0.0, 0.8, 0.0));

		switch (modelState) {
			case ModelState::Default:
				meshMarkers_default->Draw();
				break;
			case ModelState::Scaled:
				meshMarkers_scaled->Draw();
				break;
			case ModelState::Transformed:
				meshMarkers_transformedPerFrame.at(frameID)->Draw();
				break;
			case ModelState::ScapeSpace:
				// no marker data can be optained in scape space...
				break;
		}

	}
	shader->Unuse();
	glPointSize(DEFAULT_POINT_SIZE);
}


void ModelRenderer::Draw(int frameID) {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("model");

	shader->Use();
	{
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);
		shader->Uniform("cameraPosition", CameraManager::GetInstance().GetActive()->GetPosition());
		shader->Uniform("showNormals", Settings::GetInstance().showNormals);
		shader->Uniform("showBodyparts", Settings::GetInstance().showBodyParts);
		shader->Uniform("showMuscleActivity", Settings::GetInstance().showMuscles);

		switch (modelState) {
			case ModelState::Default:
				mesh_default->Draw();
				bounds = mesh_default->bounds;
				break;
			case ModelState::Scaled:
				mesh_scaled->Draw();
				bounds = mesh_scaled->bounds;
				break;
			case ModelState::Transformed:
				mesh_transformedPerFrame.at(frameID)->Draw();
				bounds = mesh_transformedPerFrame.at(frameID)->bounds;
				break;
			case ModelState::ScapeSpace:
				mesh_scapeSpacePerFrame.at(frameID)->Draw();
				bounds = mesh_scapeSpacePerFrame.at(frameID)->bounds;
				break;
		}
	}
	shader->Unuse();
}


void ModelRenderer::SetModelState(ModelState s) {
	modelState = s;

	switch (modelState) {
		case ModelState::Default:
			bounds = mesh_default->bounds;
			break;
		case ModelState::Scaled:
			bounds = mesh_scaled->bounds;
			break;
		case ModelState::Transformed:
			bounds = mesh_transformedPerFrame.at(0)->bounds;
			break;
		case ModelState::ScapeSpace:
			bounds = mesh_scapeSpacePerFrame.at(0)->bounds;
			break;
	}
}


//
// OsimRenderer
//
OsimRenderer::OsimRenderer() {
	// default
	{
		//
		// bones - default
		//
		std::vector<std::pair<glm::vec3, glm::vec3>> bones;
		for (auto& j : Model::GetInstance().osimJoints_default) {
			const Data::Joint& joint = j.second;
			bones.push_back(std::make_pair(joint.globalTransformChild[3], joint.globalTransformParent[3]));
		}
		mesh_bonesDefault = Mesh::CreateLines(bones);

		//
		// joints - default
		//
		std::vector<glm::vec3> joints;
		for (auto& b : Model::GetInstance().osimBodies_default) {
			const Data::Body& body = b.second;
			joints.push_back(body.globalTransform[3]);
		}
		mesh_jointsDefault = Mesh::CreatePoints(joints);

		//
		// markers - default
		//
		std::vector<glm::vec3> markers;
		for (auto& m : Model::GetInstance().osimMarkers_default) {
			const Data::Marker& marker = m.second;
			markers.push_back(marker.globalPosition);
		}
		mesh_markersDefault = Mesh::CreatePoints(markers);

		//
		// muscles - default
		//
		for (auto& m : Model::GetInstance().osimMuscles_default) {
			const Data::Muscle& muscle = m.second;
			if (ONLY_RENDER_SINGLE_MUSCLE != "" && muscle.name != ONLY_RENDER_SINGLE_MUSCLE) continue; // debugging
			mesh_musclesDefault.push_back(Mesh::CreateLineStrip(muscle.lineSet));
		}
	}

	// scapePose
	{
		//
		// bones - scapePose
		//
		std::vector<std::pair<glm::vec3, glm::vec3>> bones;
		for (auto& j : Model::GetInstance().osimJoints_scapePose) {
			const Data::Joint& joint = j.second;
			bones.push_back(std::make_pair(joint.globalTransformChild[3], joint.globalTransformParent[3]));
		}
		mesh_bonesScapePose = Mesh::CreateLines(bones);

		//
		// joints - scapePose
		//
		std::vector<glm::vec3> joints;
		for (auto& b : Model::GetInstance().osimBodies_scapePose) {
			const Data::Body& body = b.second;
			joints.push_back(body.globalTransform[3]);
		}
		mesh_jointsScapePose = Mesh::CreatePoints(joints);

		//
		// markers - scapePose
		//
		std::vector<glm::vec3> markers;
		for (auto& m : Model::GetInstance().osimMarkers_scapePose) {
			const Data::Marker& marker = m.second;
			markers.push_back(marker.globalPosition);
		}
		mesh_markersScapePose = Mesh::CreatePoints(markers);

		//
		// muscles - scapePose
		//
		for (auto& m : Model::GetInstance().osimMuscles_scapePose) {
			const Data::Muscle& muscle = m.second;
			if (ONLY_RENDER_SINGLE_MUSCLE != "" && muscle.name != ONLY_RENDER_SINGLE_MUSCLE) continue; // debugging
			mesh_musclesScapePose.push_back(Mesh::CreateLineStrip(muscle.lineSet));
		}
	}

	// animation
	{
		//
		// bones - animation
		//
		mesh_bonesPerFrame.resize(Model::GetInstance().osimJoints_perFrame.size());
		for (int frameID = 0; frameID < mesh_bonesPerFrame.size(); frameID++) {
			std::vector<std::pair<glm::vec3, glm::vec3>> bones;
			for (auto& j : Model::GetInstance().osimJoints_perFrame.at(frameID)) {
				const Data::Joint& joint = j.second;
				bones.push_back(std::make_pair(joint.globalTransformChild[3], joint.globalTransformParent[3]));
			}
			mesh_bonesPerFrame.at(frameID) = Mesh::CreateLines(bones);
		}

		//
		// joints - animation
		//
		mesh_jointsPerFrame.resize(Model::GetInstance().osimBodies_perFrame.size());
		for (int frameID = 0; frameID < mesh_jointsPerFrame.size(); frameID++) {
			std::vector<glm::vec3> jointCenters;
			for (auto& b : Model::GetInstance().osimBodies_perFrame.at(frameID)) {
				const Data::Body& body = b.second;
				jointCenters.push_back(body.globalTransform[3]);
			}
			mesh_jointsPerFrame.at(frameID) = Mesh::CreatePoints(jointCenters);
		}

		//
		// markers - animation
		//
		mesh_markersPerFrame.resize(Model::GetInstance().osimMarkers_perFrame.size());
		for (int frameID = 0; frameID < mesh_markersPerFrame.size(); frameID++) {
			std::vector<glm::vec3> markers;
			for (auto& m : Model::GetInstance().osimMarkers_perFrame.at(frameID)) {
				const Data::Marker& marker = m.second;
				markers.push_back(marker.globalPosition);
			}
			mesh_markersPerFrame.at(frameID) = Mesh::CreatePoints(markers);
		}

		//
		// muscles
		//
		mesh_musclesPerFrame.resize(Model::GetInstance().osimMuscles_perFrame.size());
		if (Settings::GetInstance().visualizeMuscleActivity) {
			muscleActivationColors.resize(Model::GetInstance().osimMuscles_perFrame.size());
		}
		for (int frameID = 0; frameID < mesh_musclesPerFrame.size(); frameID++) {
			for (auto& m : Model::GetInstance().osimMuscles_perFrame.at(frameID)) {
				const Data::Muscle& muscle = m.second;
				if (ONLY_RENDER_SINGLE_MUSCLE != "" && muscle.name != ONLY_RENDER_SINGLE_MUSCLE) continue; // debugging

				mesh_musclesPerFrame.at(frameID).push_back(Mesh::CreateLineStrip(muscle.lineSet));

				if (Settings::GetInstance().visualizeMuscleActivity) {
					muscleActivationColors.at(frameID).push_back(glm::vec3(muscle.activation, 0.0, 0.0));
				}
			}
		}
	}
}

OsimRenderer::~OsimRenderer() {
	SAFE_DELETE(mesh_bonesDefault);
	SAFE_DELETE(mesh_jointsDefault);
	SAFE_DELETE(mesh_markersDefault);
	for (auto& m : mesh_musclesDefault) {
		SAFE_DELETE(m);
	}

	SAFE_DELETE(mesh_bonesScapePose);
	SAFE_DELETE(mesh_jointsScapePose);
	SAFE_DELETE(mesh_markersScapePose);
	for (auto& m : mesh_musclesScapePose) {
		SAFE_DELETE(m);
	}

	for (auto& m : mesh_bonesPerFrame) {
		SAFE_DELETE(m);
	}
	for (auto& m : mesh_jointsPerFrame) {
		SAFE_DELETE(m);
	}
	for (auto& m : mesh_markersPerFrame) {
		SAFE_DELETE(m);
	}
	for (auto& meshes : mesh_musclesPerFrame) {
		for (auto& m : meshes) {
			SAFE_DELETE(m);
		}
	}
}


void OsimRenderer::Draw(int frameID) {
	switch (osimState) {
		case OsimState::Default:
			Draw_Default();
			break;
		case OsimState::Scape:
			Draw_Scape();
			break;
		case OsimState::Animation:
			Draw_Animation(frameID);
			break;
	}
}

void OsimRenderer::Draw_Default() {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("color");

	// draw OSIM data
	shader->Use();
	{
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);

		if (Settings::GetInstance().showBonesOsim) {
			// bones
			glLineWidth(2.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.7));
			mesh_bonesDefault->Draw();

			// joint centers
			glPointSize(5.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.0));
			mesh_jointsDefault->Draw();
		}

		if (Settings::GetInstance().showMarkersOsim) {
			// markers
			glPointSize(6.0f);
			shader->Uniform("color", glm::vec3(0.8, 0.0, 0.0));
			mesh_markersDefault->Draw();
		}

		if (Settings::GetInstance().showMusclesOsim) {
			// muscles
			glLineWidth(3.0f);
			for (int muscleID = 0; muscleID < mesh_musclesDefault.size(); muscleID++) {
				shader->Uniform("color", glm::vec3(0.0));
				mesh_musclesDefault.at(muscleID)->Draw();
			}
		}
	}
	shader->Unuse();
	glLineWidth(DEFAULT_LINE_WIDTH);
	glPointSize(DEFAULT_POINT_SIZE);
}


void OsimRenderer::Draw_Scape() {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("color");

	// draw OSIM data
	shader->Use();
	{
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);

		if (Settings::GetInstance().showBonesOsim) {
			// bones
			glLineWidth(2.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.7));
			mesh_bonesScapePose->Draw();

			// joint centers
			glPointSize(5.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.0));
			mesh_jointsScapePose->Draw();
		}

		if (Settings::GetInstance().showMarkersOsim) {
			// markers
			glPointSize(6.0f);
			shader->Uniform("color", glm::vec3(0.8, 0.0, 0.0));
			mesh_markersScapePose->Draw();
		}

		if (Settings::GetInstance().showMusclesOsim) {
			// muscles
			glLineWidth(3.0f);
			for (int muscleID = 0; muscleID < mesh_musclesScapePose.size(); muscleID++) {
				shader->Uniform("color", glm::vec3(0.0));
				mesh_musclesScapePose.at(muscleID)->Draw();
			}
		}
	}
	shader->Unuse();
	glLineWidth(DEFAULT_LINE_WIDTH);
	glPointSize(DEFAULT_POINT_SIZE);
}

void OsimRenderer::Draw_Animation(int frameID) {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("color");

	// draw OSIM data
	shader->Use();
	{
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);

		if (Settings::GetInstance().showBonesOsim) {
			// bones
			glLineWidth(2.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.7));
			mesh_bonesPerFrame.at(frameID)->Draw();

			// joint centers
			glPointSize(5.0f);
			shader->Uniform("color", glm::vec3(0.0, 0.7, 0.0));
			mesh_jointsPerFrame.at(frameID)->Draw();
		}

		if (Settings::GetInstance().showMarkersOsim) {
			// markers
			glPointSize(6.0f);
			shader->Uniform("color", glm::vec3(0.8, 0.0, 0.0));
			mesh_markersPerFrame.at(frameID)->Draw();
		}

		if (Settings::GetInstance().showMusclesOsim) {
			// muscles
			glLineWidth(3.0f);
			for (int muscleID = 0; muscleID < mesh_musclesPerFrame.at(frameID).size(); muscleID++) {
				if (Settings::GetInstance().visualizeMuscleActivity) {
					shader->Uniform("color", muscleActivationColors.at(frameID).at(muscleID));
				}
				mesh_musclesPerFrame.at(frameID).at(muscleID)->Draw();
			}
		}
	}
	shader->Unuse();
	glLineWidth(DEFAULT_LINE_WIDTH);
	glPointSize(DEFAULT_POINT_SIZE);
}

