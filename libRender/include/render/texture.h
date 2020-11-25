#pragma once

#include "global.h"

//
// Texture
//
class Texture {
private:
	std::string name;
	uint width, height;
	uint numMipMapLevels;

public:
	Texture(const std::string& name, GLuint target);
	~Texture();

	GLuint target;
	GLuint textureID;

	void SetTexParameter(GLenum param, GLint value);
	void SetTexImageMultisample(GLsizei samples, GLenum format, uint width, uint height, bool fixedsamplelocations);
	void SetTexImage(GLint level, GLint internalFormat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid* data);

	void Bind(int unit);
	void Bind();
	void Unbind();
};
