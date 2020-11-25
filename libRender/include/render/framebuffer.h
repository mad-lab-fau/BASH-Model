#pragma once

#include "global.h"
#include "texture.h"

#include <GL/GLU.h>
#include <GL/glew.h>


//
// Renderbuffer
//
class Renderbuffer {
private:
	std::string name;
	uint width, height;

public:
	Renderbuffer(const std::string& name, uint width, uint height);
	~Renderbuffer();

	GLuint rboID;

	void SetStorage(GLenum format);
	void SetStorageMultisample(GLenum format, GLsizei samples);

	void Bind();
};

//
// Framebuffer
//
class Framebuffer {
private:
	std::string name;
	uint width, height;

	GLint maxNumberOfColorAttachments;
	GLuint fboID;
	std::vector<Texture*> colorTextures;
	std::vector<GLenum> colorAttachments;

	Texture* depthTexture;

public:
	Framebuffer(const std::string&, uint width, uint height);

	void AttachTextureAsColor(Texture* t);
	void AttachTextureAsDepth(Texture* t);
	void AttachRenderbufferAs(Renderbuffer* r, GLenum as);

	void Bind();
	void Bind(int i);
	static void Unbind();
};
