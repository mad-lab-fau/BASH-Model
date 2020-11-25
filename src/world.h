#pragma once

#include "util.h"
#include "settings.h"
#include "model.h"
#include "renderer.h"

#include <render/camera.h>
#include <render/game.h>
#include <render/graphics.h>


//
// World
//
class World : public Renderer {
private:
	int windowWidth, windowHeight;

	int frameID = 0;
	bool incrementFrame = false;
	bool decrementFrame = false;

	// camera
	PerspectiveCamera* camera;
	OrbitCameraControl* orbitCameraControl;
	sf::Vector2i oldMousePos;
	bool camRotate = false;
	bool camPan = false;

	// render objects
	BoundingBoxRenderer* boundingBox;
	ModelRenderer* modelRenderer;
	OsimRenderer* osimRenderer;
	CoordinateSystemRenderer* originAxes;

	// Lights
	std::vector<PointLight> pointLights;
	ShaderBuffer pointLightBuffer;

public:
	World(int windowWidth, int windowHeight);
	~World();

	void Update(sf::Time delta);
	void Render(float alpha, sf::Time delta);

	void Resize(int width, int height);

	void ProcessInput();
	void ProcessMouseInput(sf::Vector2i mousePos, int mouseButton, int mouseWheel);
	void DrawFloor();

	Ray GetRayFromMousePos(glm::vec2 mousePos);
};
