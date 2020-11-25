#pragma once

#include "global.h"

//
// Bounds
//
template <typename T> // glm::vec2 | glm::vec3
class Bounds {
private:
	T bMin;
	T bMax;
	bool isValid = false;

public:
	Bounds() : bMin(T(std::numeric_limits<float>::infinity())), bMax(T(-std::numeric_limits<float>::infinity())) {};
	Bounds(const std::vector<T>& points) {
		if (points.size() > 0) {
			bMin = points[0];
			bMax = points[0];
			ExpandBy(points);
			isValid = true;
		}
	}
	T GetMin() const {
		if (isValid) {
			return bMin;
		} else {
			PRINT_ERR("GetMin() Bounds are not valid.");
			return T(0.0);
		}
	}
	T GetMax() const {
		if (isValid) {
			return bMax;
		} else {
			PRINT_ERR("GetMax() Bounds are not valid.");
			return T(0.0);
		}
	}
	T GetCenter() const {
		if (isValid) {
			return (bMax - bMin) / T(2.0) + bMin;
		} else {
			PRINT_ERR("GetCenter() Bounds are not valid.");
			return T(0.0);
		}
	}
	T GetSize() const {
		if (isValid) {
			return bMax - bMin;
		} else {
			PRINT_ERR("GetSize() Bounds are not valid.");
			return T(0.0);
		}
	}
	void ExpandBy(const T& position) {
		bMin = glm::min(bMin, position);
		bMax = glm::max(bMax, position);
		isValid = true;
	}
	void ExpandBy(const std::vector<T>& positions) {
		for (auto& p : positions) {
			ExpandBy(p);
		}
	}
	void ExpandBy(const T& min, const T& max) {
		bMin = glm::min(bMin, min);
		bMax = glm::max(bMax, max);
		isValid = true;
	}
	void ExpandBy(const Bounds<glm::vec3>& other) {
		ExpandBy(other.bMin, other.bMax);
		isValid = true;
	}
};


//
// VertexBuffer
//
class VertexBuffer {
private:
	uint index, size, dimension;
	std::string name;
	GLuint id;
	GLenum contentType, usageHint;

	void SetVertexAttrib();
	int glTypeSize(GLenum type);
public:
	VertexBuffer(uint index, const std::string& name, uint size, uint dim, GLenum contentType, GLenum usageHint);
	VertexBuffer(uint index, const std::string& name, uint size, uint dim, GLenum contentType, GLenum usageHint, const void* data);

	void ChangeData(const void* data);

	std::string GetName() { return name; }
};

//
// IndexBuffer
//
class IndexBuffer {
private:
	GLuint id;
	GLenum usageHint;

public:
	uint size;
	IndexBuffer(uint size, GLenum usageHint, const uint* data);
};

//
// Mesh
//
class Mesh {
private:
	GLuint vao;
	std::string name;
	bool bound;
	uint currentBuffer, numVertices;
	std::vector<VertexBuffer> vertexBuffers;
	IndexBuffer* indexBuffer;

	void CheckBind();
public:
	Mesh(const std::string& name);
	~Mesh();

	GLenum primitiveType;
	Bounds<glm::vec3> bounds;

	void setBound(bool b);
	void BindVAO();

	void AddBuffer(const std::string&, uint size, uint element_dim, GLenum contentType, GLenum usageHint, const void* data);
	void AddIndexBuffer(uint size, GLenum usageHint, const uint* data);

	// convinient functions
	void AddColors(const std::vector<glm::vec3>& colors);
	void AddNormals(const std::vector<glm::vec3>& normals);

	void Draw();
	void DrawAs(GLenum primitiveType);
	void Draw(uint numIndices, uint baseIndex, uint baseVertex);
	void DrawAs(GLenum primitiveType, uint numIndices, uint baseIndex, uint baseVertex);

	static Mesh* CreatePoint(const glm::vec3& vertex);
	static Mesh* CreatePoints(const std::vector<glm::vec3>& vertices);
	static Mesh* CreateLine(const Line& line);
	static Mesh* CreateLine(const glm::vec3& a, const glm::vec3& b);
	static Mesh* CreateLineStrip(const std::vector<glm::vec3>& pointSet);
	static Mesh* CreateLines(const std::vector<Line>& lines);
	static Mesh* CreateLines(const std::vector<std::pair<glm::vec3, glm::vec3>>& lines);
	static Mesh* CreateCircle(float r, uint segments);
	static Mesh* CreateRing(float r1, float r2, uint segments);
	static Mesh* CreateCirclePlane(float r, uint segments);
	static Mesh* CreatePlane(float width, float height);
	static Mesh* CreatePlane(float width, float height, uint segmentsX, uint segmentsY);
	static Mesh* CreateScreenQuad();
	static Mesh* CreateIcosphere(float r, uint subdiv);
	static Mesh* CreateSphere(float radius, uint rings, uint sectors);
	static Mesh* CreateCubeBox();
};

//
// BindManager
//
class BindManager {
private:
	Mesh* boundMesh = NULL;

	BindManager() {};
	BindManager(const BindManager&) = delete;
	BindManager& operator=(const BindManager&) = delete;
public:
	static BindManager& GetInstance() {
		static BindManager instance;
		return instance;
	}

	void Bind(Mesh* mesh);
	void Unbind();
};

