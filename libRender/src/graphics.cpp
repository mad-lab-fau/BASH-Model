#include "graphics.h"

#include "camera.h"


//
// BoundingBoxRenderer
//
BoundingBoxRenderer::BoundingBoxRenderer() {
	mesh = Mesh::CreateCubeBox();
}

BoundingBoxRenderer::~BoundingBoxRenderer() {
	SAFE_DELETE(mesh);
}

void BoundingBoxRenderer::Draw(const Bounds<glm::vec3>& bounds, const glm::mat4& M) {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().CreateColorShader();

	shader->Use();
	{
		glm::mat4 transform = glm::translate(bounds.GetCenter()) * glm::scale(bounds.GetSize());

		shader->Uniform("MVP", P * V * M * transform);
		shader->Uniform("color", color);

		mesh->Draw();
	}
	shader->Unuse();
}


//
// CoordinateSystemRenderer
//
CoordinateSystemRenderer::CoordinateSystemRenderer() {
	glm::vec3 C = glm::vec3(0, 0, 0);
	glm::vec3 X = glm::vec3(1, 0, 0);
	glm::vec3 Y = glm::vec3(0, 1, 0);
	glm::vec3 Z = glm::vec3(0, 0, 1);

	mesh_origin = Mesh::CreatePoint(C);
	mesh_axisX = Mesh::CreateLine(C, X);
	mesh_axisY = Mesh::CreateLine(C, Y);
	mesh_axisZ = Mesh::CreateLine(C, Z);
}

CoordinateSystemRenderer::~CoordinateSystemRenderer() {
	SAFE_DELETE(mesh_origin);
	SAFE_DELETE(mesh_axisX);
	SAFE_DELETE(mesh_axisY);
	SAFE_DELETE(mesh_axisZ);
}

void CoordinateSystemRenderer::Draw() {
	glm::mat4 V = CameraManager::GetInstance().GetActive()->GetV();
	glm::mat4 P = CameraManager::GetInstance().GetActive()->GetP();
	Shader* shader = ShaderManager::GetInstance().CreateColorShader();

	shader->Use();
	{
		glm::mat4 M = glm::scale(Transform, glm::vec3(lengthScale));

		shader->Uniform("MVP", P * V * M);

		shader->Uniform("color", color_origin);
		mesh_origin->Draw();

		shader->Uniform("color", color_axisX);
		mesh_axisX->Draw();

		shader->Uniform("color", color_axisY);
		mesh_axisY->Draw();

		shader->Uniform("color", color_axisZ);
		mesh_axisZ->Draw();
	}
	shader->Unuse();
}


//
// MouseCircle
//
void MouseCircle::AssignWindow(sf::RenderWindow* w) {
	window = w;
	hasWindow = true;
}

bool MouseCircle::GetVisible() {
	return visible;
}
void MouseCircle::SetVisible(bool v) {
	visible = v;
	if (hasWindow) {
		if (visible) {
			window->setMouseCursorVisible(false);
		} else {
			window->setMouseCursorVisible(true);
		}
	} else {
		PRINT_WAR("MouseCircle has no window assigned.")
	}
}

float MouseCircle::GetRadius() {
	return shape.getRadius();
}
void MouseCircle::SetRadius(float r) {
	shape.setRadius(r);
	shape.setOrigin(r, r);
}

void MouseCircle::SetPosition(float x, float y) {
	shape.setPosition(x, y);
}

void MouseCircle::Init(sf::Color color, float outlineThickness) {
	shape.setFillColor(sf::Color::Transparent);
	shape.setOutlineColor(color);
	shape.setOutlineThickness(outlineThickness);
	shape.setPosition(0.0f, 0.0f);
	shape.setRadius(0.0f);
}

void MouseCircle::Draw() {
	if (visible) {
		if (hasWindow) {
			window->setMouseCursorVisible(false);
			window->draw(shape);
		}
	}
}