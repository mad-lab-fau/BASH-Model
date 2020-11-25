#pragma once

//
// global includes for header files
//
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>


//
// Macros
//
#define PRINT(X) std::cout << X << std::endl;
#define PRINT_VAR(X) std::cout << (#X) << ": " << (X) << std::endl;
#define PRINT_ERR(X) std::cerr << "Error: " << (X) << " @ line " <<  __LINE__ << " in file: " << __FILE__ << std::endl; throw;
#define PRINT_WAR(X) std::cout << "Warning: " << (X) << " @ line " << __LINE__ << " in file: " << __FILE__ << std::endl;

#define SAFE_DELETE(p) if (p) { delete p; p = nullptr; }

//
// Window
//
namespace WindowMode {
	enum { desktop = 0, fullscreen = 1, window = 2 };
};

//
// Types
//
typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned long ulong;

struct Line {
	glm::vec3 A, B;
};
struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;
};
struct PointLight {
	glm::vec3 position;
	float shininess;
	glm::vec4 ambient;
	glm::vec4 diffuse;
	glm::vec4 specular;
};


//
// Functions
//
inline float GetDistanceRayToPoint(const Ray& ray, const glm::vec3& point) {
	return glm::length(glm::cross(ray.direction, point - ray.origin));
}

template <typename T>
inline T clamp(const T& n, const T& lower, const T& upper) {
	return glm::max(lower, glm::min(n, upper));
}

inline bool CheckOpenGLError(const char* file, int line) {
	GLenum glErr = glGetError();
	if (glErr != GL_NO_ERROR) {
		printf("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
	}
	return glErr == 0;
}
#define CHECK_OPENGL_ERROR() CheckOpenGLError(__FILE__, __LINE__);


//
// IO
//
inline std::string ReadFromFile(const std::string& filename) {
	std::string result;
	std::ifstream stream(filename, std::ios::in);

	if (stream.is_open()) {
		std::string line = "";
		while (std::getline(stream, line)) {
			result += line + "\n";
		}
		stream.close();
	} else {
		PRINT_ERR("Unable to read file " + filename);
	}
	return result;
}
