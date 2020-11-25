#pragma once

#include "global.h"

#include "shader.h"
#include "mesh.h"
#include <SFML/Graphics.hpp>

//
// BoundingBoxRenderer
//
class BoundingBoxRenderer {
private:
	Mesh* mesh;

public:
	BoundingBoxRenderer();
	~BoundingBoxRenderer();

	glm::vec3 color = glm::vec3(1.0, 0.0, 0.0);

	void Draw(const Bounds<glm::vec3>& bounds, const glm::mat4& M = glm::identity<glm::mat4>());
};

//
// CoordinateSystemRenderer
//
class CoordinateSystemRenderer {
private:
	Mesh* mesh_origin;
	Mesh* mesh_axisX;
	Mesh* mesh_axisY;
	Mesh* mesh_axisZ;

public:
	CoordinateSystemRenderer();
	~CoordinateSystemRenderer();

	float lengthScale = 0.2f;

	glm::mat4 Transform = glm::identity<glm::mat4>();
	glm::vec3 color_origin = glm::vec3(1.0, 1.0, 1.0);
	glm::vec3 color_axisX = glm::vec3(1.0, 0.0, 0.0);
	glm::vec3 color_axisY = glm::vec3(0.0, 1.0, 0.0);
	glm::vec3 color_axisZ = glm::vec3(0.0, 0.0, 1.0);

	void Draw();
};



//
// MouseCircle
//
class MouseCircle {
	sf::CircleShape shape;
	sf::RenderWindow* window = NULL;
	bool hasWindow = false;
	bool visible = false;

	MouseCircle() {};
	MouseCircle(const MouseCircle&) = delete;
	MouseCircle& operator=(const MouseCircle&) = delete;
public:
	static MouseCircle& GetInstance() {
		static MouseCircle instance;
		return instance;
	}

	bool GetVisible();
	void SetVisible(bool visible);
	float GetRadius();
	void SetRadius(float radius);
	void SetPosition(float x, float y);

	void AssignWindow(sf::RenderWindow* window);
	void Init(sf::Color color, float outlineThickness);
	void Draw();
};
