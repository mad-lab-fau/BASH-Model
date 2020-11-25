#include "camera.h"

//
// Camera
//
Camera::Camera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up)
	: name(name), position(pos), target(tar), upDir(up) {

	viewDir = glm::normalize(target - position);
	rightDir = glm::normalize(glm::cross(viewDir, upDir));

	Update();
}

Camera::~Camera() {}

void Camera::SetPosition(const glm::vec3& pos) {
	position = pos;
	viewDir = glm::normalize(target - position);
	rightDir = glm::normalize(glm::cross(viewDir, upDir));
}
void Camera::SetTarget(const glm::vec3& tar) {
	target = tar;
	viewDir = glm::normalize(target - position);
	rightDir = glm::cross(viewDir, upDir);
}
void Camera::SetUpDir(const glm::vec3& up) {
	upDir = up;
	rightDir = glm::cross(viewDir, upDir);
}
void Camera::SetPositionAndTarget(const glm::vec3& pos, const glm::vec3& tar) {
	position = pos;
	target = tar;
	viewDir = glm::normalize(target - position);
	rightDir = glm::cross(viewDir, upDir);
}
void Camera::Set(const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up) {
	position = pos;
	target = tar;
	upDir = up;
	viewDir = glm::normalize(target - position);
	rightDir = glm::cross(viewDir, upDir);
}

void Camera::MoveForward(float delta) {
	position += delta * viewDir;
	target = position + viewDir;
}
void Camera::MoveBackward(float delta) {
	position -= delta * viewDir;
	target = position + viewDir;
}
void Camera::MoveUp(float delta) {
	position += delta * upDir;
}
void Camera::MoveDown(float delta) {
	position -= delta * upDir;
}
void Camera::MoveRight(float delta) {
	position += delta * rightDir;
}
void Camera::MoveLeft(float delta) {
	position -= delta * rightDir;
}


void Camera::ScrollPosition(const glm::vec3& pos) {
	glm::vec3 diffPos = position - pos;
	position = pos;
	target -= diffPos;
}
void Camera::ScrollTarget(const glm::vec3& tar) {
	glm::vec3 diffTar = target - tar;
	target = tar;
	position -= diffTar;
}

void Camera::ScrollUp(float delta) {
	position += delta * upDir;
	target = position + viewDir;
}
void Camera::ScrollDown(float delta) {
	position -= delta * upDir;
	target = position + viewDir;
}
void Camera::ScrollRight(float delta) {
	position += delta * rightDir;
	target = position + viewDir;
}
void Camera::ScrollLeft(float delta) {
	position -= delta * rightDir;
	target = position + viewDir;
}


void Camera::Yaw(float angle) {
	viewDir = glm::normalize(glm::rotate(viewDir, angle, upDir));
	target = position + viewDir;
	rightDir = glm::cross(viewDir, upDir);
}


void Camera::Pitch(float angle) {
	viewDir = glm::normalize(glm::rotate(viewDir, angle, rightDir));
	target = position + viewDir;
	upDir = glm::cross(rightDir, viewDir);
}


void Camera::Roll(float angle) {
	rightDir = glm::normalize(glm::rotate(rightDir, angle, viewDir));
	upDir = glm::cross(rightDir, viewDir);
}


void Camera::Update() {
	V = glm::lookAt(position, target, upDir);
}

void Camera::Use() {
	if (!active) {
		CameraManager::GetInstance().UseCamera(this);
	}
}


//
// PerspectiveCamera
//
PerspectiveCamera::PerspectiveCamera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up, float fovy, int w, int h, float n, float f) :
	Camera(name, pos, tar, up), fovy(fovy), width(w), height(h), frustum_near(n), frustum_far(f) {

	glViewport(0, 0, width, height);
	active = false;
	Update();
	CameraManager::GetInstance().AddCamera(this);
}

void PerspectiveCamera::Change(float fo, int w, int h, float n, float f) {
	fovy = fo;
	width = w;
	height = h;
	frustum_near = n;
	frustum_far = f;
}

void PerspectiveCamera::Resize(int w, int h) {
	width = w;
	height = h;
	glViewport(0, 0, width, height);
	Update();
}

void PerspectiveCamera::Update() {
	Camera::Update();

	P = glm::perspective(fovy, (width * 1.0f) / height, frustum_near, frustum_far);
}


//
// OrthographicCamera
//
OrthographicCamera::OrthographicCamera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up, float left, float right, float bottom, float top, float n, float f) :
	Camera(name, pos, tar, up), left(left), right(right), bottom(bottom), top(top), frustum_near(n), frustum_far(f) {

	active = false;
	Update();
	CameraManager::GetInstance().AddCamera(this);
}

void OrthographicCamera::Update() {
	P = glm::ortho(left, right, bottom, top, frustum_near, frustum_far);
}

//
// Orbit Camera Control
//
OrbitCameraControl::OrbitCameraControl(Camera* c) : camera(c) {}

void OrbitCameraControl::SetDefault() {
	orbitCenter = glm::vec3(0.0);
	distance = 1.0;
	phi = glm::half_pi<float>();
	theta = 0.0;
}

void OrbitCameraControl::RotateToFront() {
	theta = 0.0;
	phi = glm::half_pi<float>();
	Update();
}
void OrbitCameraControl::RotateToSide() {
	theta = glm::half_pi<float>();
	phi = glm::half_pi<float>();
	Update();
}
void OrbitCameraControl::MoveOrbitCenter(const glm::vec3& c) {
	orbitCenter = c;
}
void OrbitCameraControl::Update() {
	glm::vec3 eye;
	eye.x = orbitCenter.x + distance * sin(phi) * cos(theta);
	eye.y = orbitCenter.y + distance * cos(phi) * -1;
	eye.z = orbitCenter.z + distance * sin(phi) * sin(theta);

	camera->SetPosition(eye);
	camera->SetTarget(orbitCenter);
	camera->SetUpDir(glm::vec3(0.0, 1.0, 0.0));
	camera->Update();
}


//
// Camera management system
//
Camera* CameraManager::GetCamera(const std::string& name) {
	if (cameras.find(name) == cameras.end()) {
		PRINT_ERR("Camera does not exist: " + name);
		return nullptr;
	} else {
		return cameras.at(name);
	}
}

bool CameraManager::FindCamera(const std::string& name) {
	if (cameras.find(name) == cameras.end()) {
		return false;
	} else {
		return true;
	}
}

Camera* CameraManager::GetActive() {
	return active;
}

void CameraManager::AddCamera(Camera* c) {
	if (!FindCamera(c->name)) {
		cameras[c->name] = c;
	} else {
		PRINT_ERR("Camera already exists: " + c->name);
	}
}

void CameraManager::UseCamera(Camera* c) {
	if (active) {
		active->active = false;
	}
	active = c;
	c->active = true;
}
