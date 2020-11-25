#include "framebuffer.h"


//
// Renderbuffer
//
Renderbuffer::Renderbuffer(const std::string&, uint width, uint height) :
	name(name), width(width), height(height) {

	glGenRenderbuffers(1, &rboID);
	CHECK_OPENGL_ERROR();
}

Renderbuffer::~Renderbuffer() {
	glDeleteRenderbuffers(1, &rboID);
}

void Renderbuffer::SetStorage(GLenum format) {
	glRenderbufferStorage(GL_RENDERBUFFER, format, width, height);
	CHECK_OPENGL_ERROR();
}

void Renderbuffer::SetStorageMultisample(GLenum format, GLsizei samples) {
	glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, format, width, height);
	CHECK_OPENGL_ERROR();
}

void Renderbuffer::Bind() {
	glBindRenderbuffer(GL_RENDERBUFFER, rboID);
}

//
// Framebuffer
//
Framebuffer::Framebuffer(const std::string&, uint width, uint height) :
	name(name), width(width), height(height), depthTexture(nullptr) {

	glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxNumberOfColorAttachments);

	glGenFramebuffers(1, &fboID);
	CHECK_OPENGL_ERROR();
}

void Framebuffer::AttachTextureAsColor(Texture* t) {
	int size = static_cast<int>(colorAttachments.size());
	if (size == maxNumberOfColorAttachments) {
		PRINT_ERR("Too many color attachments already");
		exit(EXIT_FAILURE);
	}
	GLenum ca = GL_COLOR_ATTACHMENT0 + size;

	colorTextures.push_back(t);
	colorAttachments.push_back(ca);

	glFramebufferTexture2D(GL_FRAMEBUFFER, ca, t->target, t->textureID, 0);
	CHECK_OPENGL_ERROR();
}

void Framebuffer::AttachTextureAsDepth(Texture* t) {
	depthTexture = t;
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, t->target, t->textureID, 0);
	CHECK_OPENGL_ERROR();
}

void Framebuffer::AttachRenderbufferAs(Renderbuffer* r, GLenum as) {
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, as, GL_RENDERBUFFER, r->rboID);
	CHECK_OPENGL_ERROR();
}

void Framebuffer::Bind() {
	glBindFramebuffer(GL_FRAMEBUFFER, fboID);
	glViewport(0, 0, width, height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
}

void Framebuffer::Bind(int i) {
	glBindFramebuffer(GL_FRAMEBUFFER, fboID);
	glViewport(0, 0, width, height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

	glDrawBuffers(static_cast<GLsizei>(1), &colorAttachments.at(i));
}

void Framebuffer::Unbind() {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
