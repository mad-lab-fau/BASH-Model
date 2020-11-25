#include "world.h"


#include <render/shader.h>
#include <render/input.h>


//
// World
//
World::World(int w, int h) : windowWidth(w), windowHeight(h) {
	// Shader
	ShaderManager::GetInstance().LoadShader("standard", "shader/standard.vert", "shader/standard.frag");
	ShaderManager::GetInstance().LoadShader("model", "shader/model.vert", "shader/model.frag");
	ShaderManager::GetInstance().LoadShader("floor", "shader/floor.vert", "shader/floor.frag");
	ShaderManager::GetInstance().LoadShader("color", "shader/color.vert", "shader/color.frag");
	ShaderManager::GetInstance().LoadShader("skinned", "shader/skinned.vert", "shader/skinned.frag");

	// Model
	modelRenderer = new ModelRenderer();

	// Osim
	osimRenderer = new OsimRenderer();

	// Bounding Box
	boundingBox = new BoundingBoxRenderer();

	// Origin & Axes
	originAxes = new CoordinateSystemRenderer();
	originAxes->Transform = glm::identity<glm::mat4>();

	// Camera
	camera = new PerspectiveCamera(CAMERA_NAME, glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0), CAMERA_FOV, windowWidth, windowHeight, 0.01f, 1000.0f);
	camera->Use();
	orbitCameraControl = new OrbitCameraControl(camera);
	orbitCameraControl->distance = CAMERA_DISTANCE;
	orbitCameraControl->orbitCenter = modelRenderer->bounds.GetCenter();
	orbitCameraControl->Update();

	// Lights
	PointLight frontLight;
	frontLight.position = POSITION_FRONT_LIGHT;
	frontLight.ambient = glm::vec4(0.05, 0.05, 0.05, 1.0);
	frontLight.diffuse = glm::vec4(0.6, 0.6, 0.6, 1.0);
	frontLight.specular = glm::vec4(0.1, 0.1, 0.1, 1.0);
	frontLight.shininess = 16;
	pointLights.push_back(frontLight);

	PointLight backLight;
	backLight.position = POSITION_BACK_LIGHT;
	backLight.ambient = glm::vec4(0.0, 0.0, 0.0, 1.0);
	backLight.diffuse = glm::vec4(0.1, 0.1, 0.1, 1.0);
	backLight.specular = glm::vec4(0.0, 0.0, 0.0, 1.0);
	backLight.shininess = 16;
	pointLights.push_back(backLight);

	pointLightBuffer.WriteData(GL_SHADER_STORAGE_BUFFER, sizeof(PointLight) * pointLights.size(), pointLights.data(), GL_STATIC_DRAW);


	// TODO: Shader buffer for materials?
}



World::~World() {
	SAFE_DELETE(camera);
	SAFE_DELETE(orbitCameraControl);
	SAFE_DELETE(boundingBox);
	SAFE_DELETE(modelRenderer);
	SAFE_DELETE(osimRenderer);
}


//
// Update
//
void World::Update(sf::Time delta) {

}


//
// Render
//
void World::Render(float alpha, sf::Time delta) {

	// camera
	camera->Use();

	// wireframe / points
	if (Settings::GetInstance().showWireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else if (Settings::GetInstance().showPoints) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	// bind light sources to shader storage buffer
	pointLightBuffer.BindBase(GL_SHADER_STORAGE_BUFFER, BIND_ID_BUFFER_POINTLIGHTS);

	// playback
	if (incrementFrame) {
		if (Settings::GetInstance().repeatPlayback) {
			frameID = (frameID + 1) % Model::GetInstance().numFrames;
		} else {
			frameID = std::min(frameID + 1, Model::GetInstance().numFrames - 1);
		}
		PRINT_VAR(frameID);
		incrementFrame = false;
	}
	if (decrementFrame) {
		frameID = std::max(0, frameID - 1);
		PRINT_VAR(frameID);
		decrementFrame = false;
	}


	// draw floor
	if (Settings::GetInstance().showFloor) {
		DrawFloor();
	}
	if (Settings::GetInstance().showOrigin) {
		// draw origin & axes
		glPointSize(6.0f);
		glLineWidth(2.0f);
		originAxes->Draw();
		glLineWidth(DEFAULT_LINE_WIDTH);
		glPointSize(DEFAULT_POINT_SIZE);
	}

	// draw osim data
	if (Settings::GetInstance().showModelOsim) {
		osimRenderer->Draw(frameID);
	}

	// markers on the model
	if (Settings::GetInstance().showMarkers) {
		modelRenderer->DrawMarkers(frameID);
	}

	if (Settings::GetInstance().showModel) {
		// render the model
		modelRenderer->Draw(frameID);
	}

	// bounding box
	if (Settings::GetInstance().showBounds) {
		glLineWidth(2.0f);
		boundingBox->Draw(modelRenderer->bounds, modelRenderer->M);
		glLineWidth(DEFAULT_LINE_WIDTH);
	}

}



// gets triggered when the windows is resized
void World::Resize(int w, int h) {
	windowWidth = w;
	windowHeight = h;
	camera->Resize(windowWidth, windowHeight);
}


//
// Input
//
void World::ProcessInput() {
	// Reload all Shader
	if (Input::GetInstance().KeyPressed(KEY_RELOAD_SHADERS)) { // "/" in EN, "#" in DE
		ShaderManager::GetInstance().Reload();
		PRINT("Shader reloaded");
	}

	// Show wireframe
	if (Input::GetInstance().KeyPressed(KEY_SHOW_WIREFRAME)) {
		Settings::GetInstance().showWireframe = !Settings::GetInstance().showWireframe;
		Settings::GetInstance().showPoints = false;
		PRINT_VAR(Settings::GetInstance().showWireframe);
	}
	// Show points
	if (Input::GetInstance().KeyPressed(KEY_SHOW_POINTS)) {
		Settings::GetInstance().showPoints = !Settings::GetInstance().showPoints;
		Settings::GetInstance().showWireframe = false;
		PRINT_VAR(Settings::GetInstance().showPoints);
	}
	// Show normals
	if (Input::GetInstance().KeyPressed(KEY_SHOW_NORMALS)) {
		Settings::GetInstance().showNormals = !Settings::GetInstance().showNormals;
		PRINT_VAR(Settings::GetInstance().showNormals);
	}

	// Show floor
	if (Input::GetInstance().KeyPressed(KEY_TOOGLE_FLOOR)) {
		Settings::GetInstance().showFloor = !Settings::GetInstance().showFloor;
		PRINT_VAR(Settings::GetInstance().showFloor);
	}
	// Show origin
	if (Input::GetInstance().KeyPressed(KEY_TOOGLE_ORIGIN)) {
		Settings::GetInstance().showOrigin = !Settings::GetInstance().showOrigin;
		PRINT_VAR(Settings::GetInstance().showOrigin);
	}
	// Show bounding box
	if (Input::GetInstance().KeyPressed(KEY_TOOGLE_BOUNDING_BOX)) {
		Settings::GetInstance().showBounds = !Settings::GetInstance().showBounds;
		PRINT_VAR(Settings::GetInstance().showBounds);
	}

	// Model State
	if (Input::GetInstance().KeyPressed(KEY_MODEL_STATE_DEFAULT)) {
		modelRenderer->SetModelState(ModelState::Default);
		PRINT("Model State: Default");
	}
	if (Input::GetInstance().KeyPressed(KEY_MODEL_STATE_SCALED)) {
		modelRenderer->SetModelState(ModelState::Scaled);
		PRINT("Model State: Scaled");
	}
	if (Input::GetInstance().KeyPressed(KEY_MODEL_STATE_TRANSFORMED)) {
		modelRenderer->SetModelState(ModelState::Transformed);
		PRINT("Model State: Transformed");
	}
	if (Input::GetInstance().KeyPressed(KEY_MODEL_STATE_SCAPESPACE)) {
		modelRenderer->SetModelState(ModelState::ScapeSpace);
		PRINT("Model State: ScapeSpace");
	}
	if (Input::GetInstance().KeyPressed(KEY_OSIM_STATE_DEFAULT)) {
		osimRenderer->SetOsimState(OsimState::Default);
		PRINT("OSIM State: Default");
	}
	if (Input::GetInstance().KeyPressed(KEY_OSIM_STATE_SCAPE)) {
		osimRenderer->SetOsimState(OsimState::Scape);
		PRINT("OSIM State: Scape");
	}
	if (Input::GetInstance().KeyPressed(KEY_OSIM_STATE_ANIMATION)) {
		osimRenderer->SetOsimState(OsimState::Animation);
		PRINT("OSIM State: Animation");
	}

	// Model: OSIM
	if (Input::GetInstance().KeyPressedContinously(KEY_MODEL_OSIM)) {
		// model
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MODEL)) {
			Settings::GetInstance().showModelOsim = !Settings::GetInstance().showModelOsim;
			if (Settings::GetInstance().showModelOsim) {
				Settings::GetInstance().showBonesOsim = true;
				Settings::GetInstance().showMarkersOsim = true;
				Settings::GetInstance().showMusclesOsim = true;
			} else {
				Settings::GetInstance().showBonesOsim = false;
				Settings::GetInstance().showMarkersOsim = false;
				Settings::GetInstance().showMusclesOsim = false;
			}
			PRINT_VAR(Settings::GetInstance().showModelOsim);
		}
		// bones
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_BONES)) {
			Settings::GetInstance().showBonesOsim = !Settings::GetInstance().showBonesOsim;
			if (Settings::GetInstance().showBonesOsim) {
				Settings::GetInstance().showModelOsim = true;
			}
			PRINT_VAR(Settings::GetInstance().showBonesOsim);
		}
		// markers
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MARKERS)) {
			Settings::GetInstance().showMarkersOsim = !Settings::GetInstance().showMarkersOsim;
			if (Settings::GetInstance().showMarkersOsim) {
				Settings::GetInstance().showModelOsim = true;
			}
			PRINT_VAR(Settings::GetInstance().showMarkersOsim);
		}
		// muscles
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MUSCLES)) {
			Settings::GetInstance().showMusclesOsim = !Settings::GetInstance().showMusclesOsim;
			if (Settings::GetInstance().showMusclesOsim) {
				Settings::GetInstance().showModelOsim = true;
			}
			PRINT_VAR(Settings::GetInstance().showMusclesOsim);
		}
	} else { // Model
		// model
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MODEL)) {
			Settings::GetInstance().showModel = !Settings::GetInstance().showModel;
			PRINT_VAR(Settings::GetInstance().showModel);
		}
		// body parts
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_BONES)) {
			Settings::GetInstance().showBodyParts = !Settings::GetInstance().showBodyParts;
			PRINT_VAR(Settings::GetInstance().showBodyParts);
			if (Settings::GetInstance().showBodyParts) {
				Settings::GetInstance().showModel = true;
			}
		}
		// markers
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MARKERS)) {
			Settings::GetInstance().showMarkers = !Settings::GetInstance().showMarkers;
			PRINT_VAR(Settings::GetInstance().showMarkers);
		}
		// muscles
		if (Input::GetInstance().KeyPressed(KEY_TOOGLE_MUSCLES)) {
			if (Settings::GetInstance().visualizeMuscleActivity) {
				Settings::GetInstance().showMuscles = !Settings::GetInstance().showMuscles;
				if (Settings::GetInstance().showMuscles) {
					Settings::GetInstance().showModel = true;
				}
				PRINT_VAR(Settings::GetInstance().showMuscles);
			} else {
				PRINT_WAR("No muscle activity data is available for visualization...")
			}
		}
	}


	// playback
	if (Input::GetInstance().KeyPressedContinously(KEY_PLAYBACK_SINGLE_FRAME)) {
		// frame by frame
		if (Input::GetInstance().KeyPressed(KEY_PLAYBACK_NEXT_FRAME)) {
			incrementFrame = true;
		}
		if (Input::GetInstance().KeyPressed(KEY_PLAYBACK_PREV_FRAME)) {
			decrementFrame = true;
		}
	} else {
		// continously
		if (Input::GetInstance().KeyPressedContinously(KEY_PLAYBACK_NEXT_FRAME)) {
			incrementFrame = true;
		}
		if (Input::GetInstance().KeyPressedContinously(KEY_PLAYBACK_PREV_FRAME)) {
			decrementFrame = true;
		}
	}
	if (Input::GetInstance().KeyPressed(KEY_TOGGLE_PLAYBACK_REPEAT)) {
		Settings::GetInstance().repeatPlayback = !Settings::GetInstance().repeatPlayback;
		PRINT_VAR(Settings::GetInstance().repeatPlayback);
	}

	if (Input::GetInstance().KeyPressed(KEY_PRESENTATION_MODE)) {
		Settings::GetInstance().presentationMode = !Settings::GetInstance().presentationMode;
		if (Settings::GetInstance().presentationMode) {
			glClearColor(PRESENTATION_COLOR);
			Settings::GetInstance().showFloor = false;
			Settings::GetInstance().showOrigin = false;
			Settings::GetInstance().showModelOsim = false;
			Settings::GetInstance().showBonesOsim = false;
			Settings::GetInstance().showMarkersOsim = false;
			Settings::GetInstance().showMusclesOsim = false;
			Settings::GetInstance().showBodyParts = false;
			Settings::GetInstance().showBounds = false;
			Settings::GetInstance().showNormals = false;
		} else {
			glClearColor(DEFAULT_CLEAR_COLOR);
		}
	}

	// defined camera states
	if (Input::GetInstance().KeyPressed(KEY_RESET_CAMERA)) {
		orbitCameraControl->SetDefault();
		orbitCameraControl->distance = CAMERA_DISTANCE;
		orbitCameraControl->orbitCenter = modelRenderer->bounds.GetCenter();
		orbitCameraControl->Update();
		PRINT("Reset camera");
	}
	if (Input::GetInstance().KeyPressed(KEY_CAMERA_FRONT)) {
		orbitCameraControl->RotateToFront();
		orbitCameraControl->orbitCenter = modelRenderer->bounds.GetCenter();
		orbitCameraControl->Update();
		PRINT("Camera to front");
	}
	if (Input::GetInstance().KeyPressed(KEY_CAMERA_SIDE)) {
		orbitCameraControl->RotateToSide();
		orbitCameraControl->orbitCenter = modelRenderer->bounds.GetCenter();
		orbitCameraControl->Update();
		PRINT("Camera to side");
	}
	if (Input::GetInstance().KeyPressed(KEY_CAMERA_DEFINED)) {
		orbitCameraControl->theta = 0.6f;
		orbitCameraControl->phi = 1.8f;
		orbitCameraControl->distance = CAMERA_DISTANCE;
		orbitCameraControl->orbitCenter = modelRenderer->bounds.GetCenter();
		orbitCameraControl->Update();
		PRINT("Camera to defined");
	}

}


void World::ProcessMouseInput(sf::Vector2i mousePos, int mouseButton, int mouseWheel) {
	// zoom camera
	if (mouseWheel != 0) {
		orbitCameraControl->distance = std::max(0.1f, orbitCameraControl->distance - (mouseWheel * CAMERA_SCROLL_SENSITIVITY));
		orbitCameraControl->Update();
	}

	// rotate camera
	if (camRotate) {
		orbitCameraControl->theta = orbitCameraControl->theta + ((mousePos.x - oldMousePos.x) * CAMERA_ROTATION_SENSITIVITY);
		orbitCameraControl->phi = clamp(orbitCameraControl->phi + ((mousePos.y - oldMousePos.y) * CAMERA_ROTATION_SENSITIVITY), glm::epsilon<float>(), glm::pi<float>() - glm::epsilon<float>());
		orbitCameraControl->Update();
	}
	if (mouseButton == MOUSE_BUTTON_ROTATE_PRESSED) {
		camRotate = true;
	}
	if (mouseButton == MOUSE_BUTTON_ROTATE_RELEASED) {
		camRotate = false;
	}
	
	// pan camera
	if (camPan) {
		glm::vec3 transX = camera->GetRightDir() * -glm::vec3((mousePos.x - oldMousePos.x) * CAMERA_PANNING_SENSITIVITY);
		glm::vec3 transY = camera->GetUpDir() * glm::vec3((mousePos.y - oldMousePos.y) * CAMERA_PANNING_SENSITIVITY);

		orbitCameraControl->orbitCenter += transX + transY;
		orbitCameraControl->Update();
	}
	if (mouseButton == MOUSE_BUTTON_PAN_PPRESSED) {
		camPan = true;
	}
	if (mouseButton == MOUSE_BUTTON_PAN_RELEASED) {
		camPan = false;
	}


	oldMousePos = mousePos;
}

void World::DrawFloor() {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().GetShader("floor");

	shader->Use();
	{
		float floorHeight = 0.0;
		glm::mat4 M = glm::identity<glm::mat4>() * glm::translate(glm::vec3(0.0, floorHeight, 0.0));
		shader->Uniform("M", M);
		shader->Uniform("V", V);
		shader->Uniform("P", P);
		shader->Uniform("cameraPosition", CameraManager::GetInstance().GetActive()->GetPosition());
		shader->Uniform("color", glm::vec3(0.7, 0.7, 0.7));

		Mesh::CreatePlane(15.f, 15.f)->Draw();
	}
	shader->Unuse();
}

Ray World::GetRayFromMousePos(glm::vec2 mousePos) {
	PerspectiveCamera* cam = (PerspectiveCamera*)CameraManager::GetInstance().GetActive();
	glm::vec4 viewport = glm::vec4(0.0, 0.0, cam->GetWidth(), cam->GetHeight());
	glm::mat4 V = cam->GetV();
	glm::mat4 P = cam->GetP();

	glm::vec3 screenPos = glm::vec3(mousePos.x, cam->GetHeight() - mousePos.y, -1.0);
	glm::vec3 nearPoint = glm::unProject(screenPos, V, P, viewport);
	screenPos = glm::vec3(mousePos.x, cam->GetHeight() - mousePos.y, 1.0);
	glm::vec3 farPoint = glm::unProject(screenPos, V, P, viewport);

	Ray ray;
	ray.origin = nearPoint;
	ray.direction = glm::normalize(farPoint - nearPoint);

	return ray;
}
