#include "mesh.h"

//
// VertexBuffer
//
VertexBuffer::VertexBuffer(uint index, const std::string& name, uint size, uint dim, GLenum contentType, GLenum usageHint) :
	VertexBuffer(index, name, size, dim, contentType, usageHint, 0) {}

VertexBuffer::VertexBuffer(uint index, const std::string& name, uint size, uint dim, GLenum contentType, GLenum usageHint, const void* data) :
	index(index), name(name), size(size),
	dimension(dim), contentType(contentType), usageHint(usageHint) {

	glGenBuffers(1, &id); // Performance when each buffer has seperate id?

	int typeSize = glTypeSize(contentType);
	uint sizeInBytes = typeSize * size * dimension;

	glBindBuffer(GL_ARRAY_BUFFER, id);
	glBufferData(GL_ARRAY_BUFFER, sizeInBytes, data, usageHint);
	SetVertexAttrib();
	glEnableVertexAttribArray(index);
}

void VertexBuffer::SetVertexAttrib() {
	if (contentType == GL_FLOAT) {
		glVertexAttribPointer(index, dimension, contentType, GL_FALSE, 0, 0);
	} else {
		glVertexAttribIPointer(index, dimension, contentType, 0, 0);
	}
}

void VertexBuffer::ChangeData(const void* data) {
	uint sizeInBytes = glTypeSize(contentType) * size * dimension;
	glBindBuffer(GL_ARRAY_BUFFER, id);
	glBufferData(GL_ARRAY_BUFFER, sizeInBytes, data, usageHint);
	SetVertexAttrib();
	glEnableVertexAttribArray(index);
}

int VertexBuffer::glTypeSize(GLenum type) {
	if (type == GL_FLOAT)
		return sizeof(GLfloat);
	else if (type == GL_INT)
		return sizeof(int);
	else if (type == GL_UNSIGNED_INT)
		return sizeof(GLuint);
	else {
		PRINT_ERR("glTypeSize");
		return 0;
	}
}

//
// IndexBuffer
//
IndexBuffer::IndexBuffer(uint size, GLenum usageHint, const uint* data) :
	usageHint(usageHint), size(size) {

	glGenBuffers(1, &id);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint) * size, data, usageHint);
}


//
// Mesh
//
Mesh::Mesh(const std::string& name) :
	name(name), currentBuffer(0), numVertices(0), indexBuffer(nullptr), bound(false), primitiveType(GL_TRIANGLES) {

	glGenVertexArrays(1, &vao);
}

Mesh::~Mesh() {
	SAFE_DELETE(indexBuffer);
}

void Mesh::CheckBind() {
	if (!bound) {
		BindManager::GetInstance().Bind(this);
	}
}

void Mesh::BindVAO() {
	glBindVertexArray(vao);
}

void Mesh::setBound(bool b) {
	bound = b;
}

void Mesh::AddBuffer(const std::string& name, uint size, uint dim, GLenum contentType, GLenum usageHint, const void* data) {
	if (size == 0) {
		PRINT_ERR("Empty vertex-attribute buffer: " + name);
	}
	if (this->numVertices == 0) {
		this->numVertices = size;
	} else {
		if (this->numVertices != size) {
			PRINT_ERR("Size doesn't match other VBs");
		}
	}

	CheckBind();
	VertexBuffer v(currentBuffer, name, size, dim, contentType, usageHint, data);
	vertexBuffers.push_back(v);
	currentBuffer++;
}

void Mesh::AddColors(const std::vector<glm::vec3>& colors) {
	AddBuffer("in_color", static_cast<uint>(colors.size()), 3, GL_FLOAT, GL_STATIC_DRAW, colors.data());
}

void Mesh::AddNormals(const std::vector<glm::vec3>& normals) {
	AddBuffer("in_normal", static_cast<uint>(normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, normals.data());
}

void Mesh::AddIndexBuffer(uint size, GLenum usageHint, const uint* data) {
	if (size == 0) {
		PRINT_ERR("Empty indexbuffer: " + name);
	}
	CheckBind();
	indexBuffer = new IndexBuffer(size, usageHint, data);
}


void Mesh::Draw() {
	DrawAs(this->primitiveType);
}

void Mesh::DrawAs(GLenum primitiveType) {
	CheckBind();
	if (indexBuffer) {
		glDrawElements(primitiveType, indexBuffer->size, GL_UNSIGNED_INT, 0);
	} else {
		glDrawArrays(primitiveType, 0, numVertices);
	}
}

void Mesh::Draw(uint numIndices, uint baseIndex, uint baseVertex) {
	DrawAs(this->primitiveType, numIndices, baseIndex, baseVertex);
}

void Mesh::DrawAs(GLenum primitiveType, uint numIndices, uint baseIndex, uint baseVertex) {
	CheckBind();
	if (indexBuffer) {
		glDrawElementsBaseVertex(primitiveType, numIndices, GL_UNSIGNED_INT, (void*)(sizeof(uint) * baseIndex), baseVertex);
	} else {
		PRINT_ERR("glDrawElementsBaseVertex called, but mesh has no indexBuffer.");
	}
}

//
// BindManager
//
void BindManager::Bind(Mesh* m) {
	if (m == boundMesh) {
		return;
	}
	if (boundMesh != 0) {
		boundMesh->setBound(false);
	}
	m->BindVAO();
	boundMesh = m;
	m->setBound(true);
}

void BindManager::Unbind() {
	if (boundMesh != 0) {
		boundMesh->setBound(false);
		boundMesh = 0;
		glBindVertexArray(0);
	}
}


//
// Mesh Functions
//
Mesh* Mesh::CreatePoint(const glm::vec3& p) {
	Mesh* m = new Mesh("Point");
	m->primitiveType = GL_POINTS;

	m->AddBuffer("in_position", 1, 3, GL_FLOAT, GL_STATIC_DRAW, &p);
	return m;
}

Mesh* Mesh::CreatePoints(const std::vector<glm::vec3>& vertices) {
	Mesh* m = new Mesh("Point");
	m->primitiveType = GL_POINTS;

	m->AddBuffer("in_position", vertices.size(), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}

Mesh* Mesh::CreateLine(const Line& line) {
	return CreateLine(line.A, line.B);
}
Mesh* Mesh::CreateLine(const glm::vec3& a, const glm::vec3& b) {
	Mesh* m = new Mesh("Line");
	m->primitiveType = GL_LINES;

	std::vector<glm::vec3> vertices;
	vertices.push_back(a);
	vertices.push_back(b);

	m->AddBuffer("in_position", vertices.size(), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}

Mesh* Mesh::CreateLineStrip(const std::vector<glm::vec3>& pointSet) {
	Mesh* m = new Mesh("LineStrip");
	m->primitiveType = GL_LINE_STRIP;
	m->AddBuffer("in_position", pointSet.size(), 3, GL_FLOAT, GL_STATIC_DRAW, pointSet.data());
	return m;
}

Mesh* Mesh::CreateLines(const std::vector<Line>& lines) {
	Mesh* m = new Mesh("Lines");
	m->primitiveType = GL_LINES;

	std::vector<glm::vec3> vertices;
	for (const auto& line : lines) {
		vertices.push_back(line.A);
		vertices.push_back(line.B);
	}
	m->AddBuffer("in_position", vertices.size(), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}
Mesh* Mesh::CreateLines(const std::vector<std::pair<glm::vec3, glm::vec3>>& lines) {
	Mesh* m = new Mesh("Lines");
	m->primitiveType = GL_LINES;

	std::vector<glm::vec3> vertices;
	for (const auto& line : lines) {
		vertices.push_back(line.first);
		vertices.push_back(line.second);
	}
	m->AddBuffer("in_position", vertices.size(), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}

Mesh* Mesh::CreateCircle(float r, uint segments) {
	Mesh* m = new Mesh("Circle");
	m->primitiveType = GL_LINE_LOOP;

	std::vector<glm::vec3> vertices;

	float step = glm::two_pi<float>() / segments;
	for (uint i = 0; i < segments; ++i) {
		vertices.push_back(glm::vec3(r * cos(i * step), 0, r * sin(i * step)));
	}

	m->AddBuffer("in_position", vertices.size(), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}

Mesh* Mesh::CreateRing(float r1, float r2, uint segments) {
	Mesh* m = new Mesh("Ring");
	m->primitiveType = GL_TRIANGLE_STRIP;

	std::vector<glm::vec3> vertices;

	float step = glm::two_pi<float>() / segments;
	for (uint i = 0; i <= segments; ++i) {
		vertices.push_back(glm::vec3(r1 * cos(i * step), 0, r1 * sin(i * step)));
		vertices.push_back(glm::vec3(r2 * cos(i * step), 0, r2 * sin(i * step)));
	}

	m->AddBuffer("in_position", (segments + 1) * 2, 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	return m;
}

Mesh* Mesh::CreateCirclePlane(float r, uint segments) {
	Mesh* m = new Mesh("CirclePlane");
	m->primitiveType = GL_TRIANGLE_FAN;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;

	float step = glm::two_pi<float>() / segments;
	for (uint i = 0; i <= segments; ++i) {
		vertices.push_back(glm::vec3(r * cos(i * step), 0, r * sin(i * step)));
		normals.push_back(glm::vec3(0.0f, -1.0f, 0.0f));
	}

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddBuffer("in_normal", static_cast<uint>(normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, normals.data());
	return m;
}

Mesh* Mesh::CreatePlane(float width, float height) {
	Mesh* m = new Mesh("Plane");
	m->primitiveType = GL_TRIANGLES;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;

	float w = width / 2;
	float h = height / 2;
	vertices = {
		{ w, 0.0f, h },
		{ -w, 0.0f, h },
		{ -w, 0.0f, -h },
		{ -w, 0.0f, -h },
		{ w, 0.0f, -h },
		{ w, 0.0f, h },
	};

	for (int i = 0; i < 6; i++) {
		normals.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	}

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddBuffer("in_normal", static_cast<uint>(normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, normals.data());

	return m;
}

Mesh* Mesh::CreatePlane(float width, float height, uint pointsX, uint pointsY) {
	Mesh* m = new Mesh("Plane");
	m->primitiveType = GL_TRIANGLES;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<uint> indices;

	float stepX = width / (pointsX - 1);
	float stepY = height / (pointsY - 1);
	for (uint y = 0; y < pointsY; ++y) {

		for (uint x = 0; x < pointsX; ++x) {
			glm::vec3 vec;
			if ((y % 2) == 0) {
				vec = glm::vec3(x * stepX, 0.0f, y * stepY);
				if (y < pointsY - 1 && x < pointsX - 1) {
					indices.push_back(y * pointsX + x);
					indices.push_back((y + 1) * pointsX + (x + 1));
					indices.push_back((y + 1) * pointsX + x);

					indices.push_back(y * pointsX + x);
					indices.push_back(y * pointsX + (x + 1));
					indices.push_back((y + 1) * pointsX + (x + 1));
				}
			} else {
				vec = glm::vec3(x * stepX - stepX / 2.0, 0.0f, y * stepY);
				if (y < pointsY - 1 && x < pointsX - 1) {
					indices.push_back(y * pointsX + x);
					indices.push_back(y * pointsX + (x + 1));
					indices.push_back((y + 1) * pointsX + x);

					indices.push_back(y * pointsX + (x + 1));
					indices.push_back((y + 1) * pointsX + (x + 1));
					indices.push_back((y + 1) * pointsX + x);
				}
			}
			vertices.push_back(vec);
			normals.push_back(glm::vec3(0.0f, 1.0f, 0.0f));


		}
	}

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddBuffer("in_normal", static_cast<uint>(normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, normals.data());
	m->AddIndexBuffer(static_cast<uint>(indices.size()), GL_STATIC_DRAW, indices.data());

	return m;
}

Mesh* Mesh::CreateScreenQuad() {
	Mesh* screenQuad = new Mesh("ScreenQuad");
	std::vector<glm::vec3> vertices = { {0.0f,0.0f,10.0f}, {1.0f,0.0f,10.0f},
										{1.0f,1.0f,10.0f}, {0.0f,1.0f,10.0f} }; // TODO: why z = 10?, quad gets clipped...
	std::vector<glm::vec2> textureCoords = { {0.0f,0.0f}, {1.0f,0.0f}, {1.0f,1.0f}, {0.0f,1.0f} };
	std::vector<uint> indices = { 0, 1, 2, 2, 3, 0 };
	screenQuad->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	screenQuad->AddBuffer("in_tc", static_cast<uint>(textureCoords.size()), 2, GL_FLOAT, GL_STATIC_DRAW, textureCoords.data());
	screenQuad->AddIndexBuffer(static_cast<uint>(indices.size()), GL_STATIC_DRAW, indices.data());

	return screenQuad;
}

uint addVertex(std::vector<glm::vec3>& vec, std::vector<glm::vec3>& normals, glm::vec3 p, float r, uint& index) {
	p = glm::normalize(p);
	normals.push_back(p);
	vec.push_back(p * r);
	return index++;
}

void addFace(std::vector<uint>& vec, uint i1, uint i2, uint i3) {
	vec.push_back(i1);
	vec.push_back(i2);
	vec.push_back(i3);
}

uint getMiddlePoint(uint p1, uint p2, std::unordered_map<uint, uint>& cache, std::vector<glm::vec3>& vertices, std::vector<glm::vec3>& normals, float r, uint& index) {
	// first check if we have it already
	bool firstIsSmaller = p1 < p2;
	uint smallerIndex = firstIsSmaller ? p1 : p2;
	uint greaterIndex = firstIsSmaller ? p2 : p1;
	uint key = (smallerIndex << 16) + greaterIndex;

	if (cache.find(key) != cache.end()) {
		return cache.at(key);
	}

	// not in cache, calculate it
	glm::vec3 point1 = vertices.at(p1);
	glm::vec3 point2 = vertices.at(p2);
	glm::vec3 middle = glm::vec3(
		(point1.x + point2.x) / 2.0f,
		(point1.y + point2.y) / 2.0f,
		(point1.z + point2.z) / 2.0f);

	// add vertex makes sure point is on unit sphere

	int i = addVertex(vertices, normals, middle, r, index); ;

	// store it, return index
	cache[key] = i;
	return i;
}

Mesh* Mesh::CreateIcosphere(float r, uint subdiv) {
	Mesh* m = new Mesh("Icosphere");

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<uint> indices;
	uint index = 0;

	float t = (1.0f + sqrt(5.0f)) / 2.0f;
	addVertex(vertices, normals, glm::vec3(-1, t, 0), r, index);
	addVertex(vertices, normals, glm::vec3(1, t, 0), r, index);
	addVertex(vertices, normals, glm::vec3(-1, -t, 0), r, index);
	addVertex(vertices, normals, glm::vec3(1, -t, 0), r, index);

	addVertex(vertices, normals, glm::vec3(0, -1, t), r, index);
	addVertex(vertices, normals, glm::vec3(0, 1, t), r, index);
	addVertex(vertices, normals, glm::vec3(0, -1, -t), r, index);
	addVertex(vertices, normals, glm::vec3(0, 1, -t), r, index);

	addVertex(vertices, normals, glm::vec3(t, 0, -1), r, index);
	addVertex(vertices, normals, glm::vec3(t, 0, 1), r, index);
	addVertex(vertices, normals, glm::vec3(-t, 0, -1), r, index);
	addVertex(vertices, normals, glm::vec3(-t, 0, 1), r, index);

	// 5 faces around point 0
	addFace(indices, 0, 11, 5);
	addFace(indices, 0, 5, 1);
	addFace(indices, 0, 1, 7);
	addFace(indices, 0, 7, 10);
	addFace(indices, 0, 10, 11);

	// 5 adjacent faces 
	addFace(indices, 1, 5, 9);
	addFace(indices, 5, 11, 4);
	addFace(indices, 11, 10, 2);
	addFace(indices, 10, 7, 6);
	addFace(indices, 7, 1, 8);

	// 5 faces around point 3
	addFace(indices, 3, 9, 4);
	addFace(indices, 3, 4, 2);
	addFace(indices, 3, 2, 6);
	addFace(indices, 3, 6, 8);
	addFace(indices, 3, 8, 9);

	// 5 adjacent faces 
	addFace(indices, 4, 9, 5);
	addFace(indices, 2, 4, 11);
	addFace(indices, 6, 2, 10);
	addFace(indices, 8, 6, 7);
	addFace(indices, 9, 8, 1);

	std::unordered_map<uint, uint> cache;

	for (uint i = 0; i < subdiv; i++) {
		uint nFaces = static_cast<uint>(indices.size());
		for (uint j = 0; j < nFaces; j += 3) {
			uint v1 = indices.at(j);
			uint v2 = indices.at(j + 1);
			uint v3 = indices.at(j + 2);

			uint a = getMiddlePoint(v1, v2, cache, vertices, normals, r, index);
			uint b = getMiddlePoint(v2, v3, cache, vertices, normals, r, index);
			uint c = getMiddlePoint(v3, v1, cache, vertices, normals, r, index);

			addFace(indices, v1, a, c);
			addFace(indices, v2, b, a);
			addFace(indices, v3, c, b);

			indices[j] = a;
			indices[j + 1] = b;
			indices[j + 2] = c;
		}
	}

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddBuffer("in_normal", static_cast<uint>(normals.size()), 3, GL_FLOAT, GL_STATIC_DRAW, normals.data());
	m->AddIndexBuffer(static_cast<uint>(indices.size()), GL_STATIC_DRAW, indices.data());

	return m;
}



Mesh* Mesh::CreateSphere(float radius, uint rings, uint sectors) {
	Mesh* m = new Mesh("Sphere");
	m->primitiveType = GL_QUADS;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> texCoords;
	std::vector<uint> indices;

	// create sphere geometry
	const float R = 1.0f / (float)(rings - 1);
	const float S = 1.0f / (float)(sectors - 1);

	for (uint r = 0; r < rings; ++r) {
		for (uint s = 0; s < sectors; ++s) {
			const float y = sin(-glm::half_pi<float>() + glm::pi<float>() * r * R);
			const float x = cos(2 * glm::pi<float>() * s * S) * sin(glm::pi<float>() * r * R);
			const float z = sin(2 * glm::pi<float>() * s * S) * sin(glm::pi<float>() * r * R);

			vertices.push_back(glm::vec3(x, y, z) * radius);
			texCoords.push_back(glm::vec2(s * S, r * R));

			if (r < rings - 1) {
				int curRow = r * sectors;
				int nextRow = (r + 1) * sectors;
				int nextS = (s + 1) % sectors;

				indices.push_back(curRow + s);
				indices.push_back(nextRow + s);
				indices.push_back(nextRow + nextS);

				indices.push_back(curRow + s);
				indices.push_back(nextRow + nextS);
				indices.push_back(curRow + nextS);
			}
		}
	}

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddBuffer("in_tc", static_cast<uint>(texCoords.size()), 2, GL_FLOAT, GL_STATIC_DRAW, texCoords.data());
	m->AddIndexBuffer(static_cast<uint>(indices.size()), GL_STATIC_DRAW, indices.data());

	return m;
}

Mesh* Mesh::CreateCubeBox() {
	Mesh* m = new Mesh("CubeBox");
	m->primitiveType = GL_LINES;

	std::vector<glm::vec3> vertices = {
		{-0.5, -0.5, -0.5},
		{0.5, -0.5, -0.5},
		{0.5, 0.5, -0.5},
		{-0.5, 0.5, -0.5},
		{-0.5, -0.5, 0.5},
		{0.5, -0.5, 0.5},
		{0.5, 0.5, 0.5},
		{-0.5, 0.5, 0.5}
	};

	std::vector<uint> indices = {
		0, 1, 1, 2, 2, 3, 3, 0, // front
		4, 5, 5, 6, 6, 7, 7, 4, // back
		0, 4, 1, 5, 2, 6, 3, 7 // rest
	};

	m->AddBuffer("in_position", static_cast<uint>(vertices.size()), 3, GL_FLOAT, GL_STATIC_DRAW, vertices.data());
	m->AddIndexBuffer(static_cast<uint>(indices.size()), GL_STATIC_DRAW, indices.data());

	return m;
}

