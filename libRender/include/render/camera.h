#pragma once

#include "global.h"

// NOTE: don't use far/near as variable names

//
// Camera
//
class Camera { // Super class for Cameras. Name, direction, movement, and camera managment system
private:
	glm::vec3 position, target, viewDir, upDir, rightDir;
	glm::mat4 V;

public:
	Camera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up);
	virtual ~Camera();

	std::string name;

	inline glm::vec3 GetPosition() { return position; }
	inline glm::vec3 GetTarget() { return target; }
	inline glm::vec3 GetViewDir() { return viewDir; }
	inline glm::vec3 GetUpDir() { return upDir; }
	inline glm::vec3 GetRightDir() { return rightDir; }

	void SetPosition(const glm::vec3& pos);
	void SetTarget(const glm::vec3& tar);
	void SetUpDir(const glm::vec3& up);
	void SetPositionAndTarget(const glm::vec3& pos, const glm::vec3& tar);
	void Set(const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up);

	void MoveForward(float delta);
	void MoveBackward(float delta);
	void MoveRight(float delta);
	void MoveLeft(float delta);
	void MoveUp(float delta);
	void MoveDown(float delta);
	void ScrollPosition(const glm::vec3& pos);
	void ScrollTarget(const glm::vec3& tar);
	void ScrollRight(float delta);
	void ScrollLeft(float delta);
	void ScrollUp(float delta);
	void ScrollDown(float delta);
	void Yaw(float angle);
	void Pitch(float angle);
	void Roll(float angle);

	virtual void Update();
	inline glm::mat4 GetV() { return V; }
	virtual glm::mat4 GetP() = 0;

	// camera management system
	bool active;
	void Use();
};


//
// PerspectiveCamera
//
class PerspectiveCamera : public Camera {
private:
	float fovy;
	int width, height;
	float frustum_near, frustum_far;
	glm::mat4 P;

public:
	PerspectiveCamera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up, float fovy, int width, int height, float frustum_near, float frustum_far);

	void Change(float fovy, int width, int height, float frustum_near, float frustum_far);
	void Resize(int width, int height);

	void Update();
	inline glm::mat4 GetP() { return P; }

	inline float GetFOV() { return fovy; }
	inline int GetWidth() { return width; }
	inline int GetHeight() { return height; }
	inline glm::vec2 GetWindowSize() { return glm::vec2(width, height); }
};

//
// OrthographicCamera
//
class OrthographicCamera : public Camera {
private:
	float left, right, bottom, top;
	float frustum_near, frustum_far;
	glm::mat4 P;
public:
	OrthographicCamera(const std::string& name, const glm::vec3& pos, const glm::vec3& tar, const glm::vec3& up, float left, float right, float bottom, float top, float frustum_near, float frustum_far);

	void Update();
	inline glm::mat4 GetP() { return P; }
};

//
// OrbitCameraControl
//
class OrbitCameraControl {
private:
	Camera* camera = NULL;

public:
	glm::vec3 orbitCenter = glm::vec3(0.0);
	float distance = 1.0;
	float phi = glm::half_pi<float>(); // yaw
	float theta = 0.0; // pitch

	OrbitCameraControl(Camera* camera);

	void SetDefault();
	void RotateToFront();
	void RotateToSide();
	void MoveOrbitCenter(const glm::vec3& c);
	void Update();
};


//
// CameraManager
//
class CameraManager {
	std::unordered_map<std::string, Camera*> cameras;
	Camera* active = NULL;

	CameraManager() {};
	CameraManager(const CameraManager&) = delete;
	CameraManager& operator=(const CameraManager&) = delete;
public:
	static CameraManager& GetInstance() {
		static CameraManager instance;
		return instance;
	}

	Camera* GetCamera(const std::string& name);
	bool FindCamera(const std::string& name);
	void AddCamera(Camera* cam);
	Camera* GetActive();
	void UseCamera(Camera* cam);
};
