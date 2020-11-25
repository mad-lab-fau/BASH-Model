#include "texture.h"


//
// Texture
//
Texture::Texture(const std::string& name, GLuint target) :
	name(name), target(target) {

	glGenTextures(1, &textureID);
}

void Texture::SetTexParameter(GLenum param, GLint value) {
	glTexParameteri(target, param, value);
}

void Texture::SetTexImageMultisample(GLsizei samples, GLenum format, uint w, uint h, bool fixedsamplelocations) {
	width = w;
	height = h;
	glTexImage2DMultisample(target, samples, format, width, height, fixedsamplelocations);
	CHECK_OPENGL_ERROR();
}

void Texture::SetTexImage(GLint level, GLint internalFormat, GLsizei w, GLsizei h, GLint border, GLenum format, GLenum type, const GLvoid* data) {
	width = w;
	height = h;
	glTexImage2D(target, level, internalFormat, width, height,
		border, format, type, data);
}

//Texture::Texture(std::string name, GLuint target) :
//		name(name), target(target), internalFormat(internalFormat), 
//		format(format), width(width), height(height) {
//
//	glGenTextures(1, &textureID);
//
//	glBindTexture(target, textureID);
//
//	SetTextureParams(params);
//
//	if (target == GL_TEXTURE_2D) {
//		glTexImage2D(target, 0, internalFormat, width, height, 0, format, type, 0);
//	} else if (target == GL_TEXTURE_2D_MULTISAMPLE) {
//		glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, internalFormat, width, height, GL_FALSE);
//	}
//	isValid = true;
//
//	glBindTexture(target, 0);
//	CHECK_OPENGL_ERROR();
//}
//
//Texture::Texture(std::string name, GLuint target, GLenum internalFormat, GLenum type, GLenum format, TextureParams params, std::string filepath) :
//		name(name), target(target), internalFormat(internalFormat), format(format) {
//
//	sf::Image image;
//	if (!image.loadFromFile(filepath)) {
//		PRINT_ERR("  Could not load texture: " + filepath);
//		return;
//	} else {
//		PRINT("  Loaded texture: " + filepath);
//		isValid = true;
//	}
//	width = image.getSize().x;
//	height = image.getSize().y;
//
//	glGenTextures(1, &textureID);
//	glBindTexture(target, textureID);
//
//	numMipMapLevels = params.mipMapping ? 1 + floor(log2(std::max(width, height))) : 1;
//	glTexStorage2D(target, numMipMapLevels, GL_RGBA8, width, height);
//	glTexSubImage2D(target, 0, 0, 0, width, height, format, type, image.getPixelsPtr());
//
//	SetTextureParams(params);
//
//	glBindTexture(target, 0);
//	CHECK_OPENGL_ERROR();
//}
//
//Texture::Texture(std::string name, GLuint target, GLenum internalFormat, GLenum type, GLenum format, TextureParams params, uint width, uint height, const void* pixels) :
//		name(name), target(target), internalFormat(internalFormat), format(format), width(width), height(height)  {
//	glGenTextures(1, &textureID);
//	glBindTexture(target, textureID);
//
//	numMipMapLevels = params.mipMapping ? 1 + floor(log2(std::max(width, height))) : 1;
//	glTexStorage2D(target, numMipMapLevels, GL_RGBA8, width, height);
//	glTexSubImage2D(target, 0, 0, 0, width, height, format, type, pixels);
//
//	SetTextureParams(params);
//
//	isValid = true;
//
//	glBindTexture(target, 0);
//	CHECK_OPENGL_ERROR();
//}

Texture::~Texture() {}

//void Texture::SetTextureParams(TextureParams params) {
//	if (params.mipMapping) {
//		glGenerateMipmap(target);
//	}
//	glTexParameteri(target, GL_TEXTURE_WRAP_S, params.wrapS);
//	glTexParameteri(target, GL_TEXTURE_WRAP_T, params.wrapT);
//	glTexParameteri(target, GL_TEXTURE_MAG_FILTER, params.mag);
//	glTexParameteri(target, GL_TEXTURE_MIN_FILTER, params.min);
//}

void Texture::Bind(int unit) {
	glActiveTexture(GL_TEXTURE0 + unit);
	glBindTexture(target, textureID);
}

void Texture::Bind() {
	glBindTexture(target, textureID);
}

void Texture::Unbind() {
	glBindTexture(target, 0);
}
