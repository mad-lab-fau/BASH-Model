#pragma once

#include "global.h"

#include <glm/gtc/type_ptr.hpp>
#define GLM_FORCE_RADIANS

//
// Shader
//
class Shader {
private:
	GLuint program;

public:
	Shader(const std::string& name);

	std::string name;

	void Use();
	void Unuse();
	virtual void Reload() = 0;

	void CompileShader(const std::string& vertexShaderCode, const std::string& fragmentShaderCode);
	void CompileShader(const std::string& vertexShaderCode, const std::string& geometryShaderCode, const std::string& fragmentShaderCode);
	static GLuint CompileAndCheck(const std::string& shaderCode, GLenum shaderType);
	static GLuint LinkAndCheck(GLuint vertexShaderID, GLuint fragmentShaderID);
	static GLuint LinkAndCheck(GLuint vertexShaderID, GLuint geometryShaderID, GLuint fragmentShaderID);

	inline int GetLoc(const std::string& uniform) { return glGetUniformLocation(program, uniform.c_str()); }

	// bool
	inline void Uniform(const std::string& uniform, bool x) { glUniform1i(GetLoc(uniform), x); }

	// float
	inline void Uniform(const std::string& uniform, float x) { glUniform1f(GetLoc(uniform), x); }
	inline void Uniform(const std::string& uniform, float x, float y) { glUniform2f(GetLoc(uniform), x, y); }
	inline void Uniform(const std::string& uniform, float x, float y, float z) { glUniform3f(GetLoc(uniform), x, y, z); }
	inline void Uniform(const std::string& uniform, float x, float y, float z, float w) { glUniform4f(GetLoc(uniform), x, y, z, w); }
	inline void Uniform(const std::string& uniform, std::vector<float> x) { glUniform1fv(GetLoc(uniform), static_cast<GLsizei>(x.size()), (float*)&x[0]); }

	// int
	inline void Uniform(const std::string& uniform, int x) { glUniform1i(GetLoc(uniform), x); }
	inline void Uniform(const std::string& uniform, int x, int y) { glUniform2i(GetLoc(uniform), x, y); }
	inline void Uniform(const std::string& uniform, int x, int y, int z) { glUniform3i(GetLoc(uniform), x, y, z); }
	inline void Uniform(const std::string& uniform, int x, int y, int z, int w) { glUniform4i(GetLoc(uniform), x, y, z, w); }
	inline void Uniform(const std::string& uniform, std::vector<int> x) { glUniform1iv(GetLoc(uniform), static_cast<GLsizei>(x.size()), (int*)&x[0]); }

	// uint
	inline void Uniform(const std::string& uniform, uint x) { glUniform1ui(GetLoc(uniform), x); }
	inline void Uniform(const std::string& uniform, uint x, uint y) { glUniform2ui(GetLoc(uniform), x, y); }
	inline void Uniform(const std::string& uniform, uint x, uint y, uint z) { glUniform3ui(GetLoc(uniform), x, y, z); }
	inline void Uniform(const std::string& uniform, uint x, uint y, uint z, uint w) { glUniform4ui(GetLoc(uniform), x, y, z, w); }

	// glm
	inline void Uniform(const std::string& uniform, const glm::vec2& vec) { Uniform(uniform, vec.x, vec.y); }
	inline void Uniform(const std::string& uniform, const glm::vec3& vec) { Uniform(uniform, vec.x, vec.y, vec.z); }
	inline void Uniform(const std::string& uniform, const glm::vec4& vec) { Uniform(uniform, vec.x, vec.y, vec.z, vec.w); }
	inline void Uniform(const std::string& uniform, const glm::mat2& mat) { glUniformMatrix2fv(GetLoc(uniform), 1, false, glm::value_ptr(mat)); }
	inline void Uniform(const std::string& uniform, const glm::mat3& mat) { glUniformMatrix3fv(GetLoc(uniform), 1, false, glm::value_ptr(mat)); }
	inline void Uniform(const std::string& uniform, const glm::mat4& mat) { glUniformMatrix4fv(GetLoc(uniform), 1, false, glm::value_ptr(mat)); }

	// vectors
	inline void Uniform(const std::string& uniform, const std::vector<glm::vec2>& vec) { glUniform2fv(GetLoc(uniform), static_cast<GLsizei>(vec.size()), (float*)&vec[0]); }
	inline void Uniform(const std::string& uniform, const std::vector<glm::vec3>& vec) { glUniform3fv(GetLoc(uniform), static_cast<GLsizei>(vec.size()), (float*)&vec[0]); }
	inline void Uniform(const std::string& uniform, const std::vector<glm::vec4>& vec) { glUniform4fv(GetLoc(uniform), static_cast<GLsizei>(vec.size()), (float*)&vec[0]); }
	inline void Uniform(const std::string& uniform, const std::vector<glm::mat2>& mat) { glUniformMatrix2fv(GetLoc(uniform), static_cast<GLsizei>(mat.size()), false, (float*)&mat[0]); }
	inline void Uniform(const std::string& uniform, const std::vector<glm::mat3>& mat) { glUniformMatrix3fv(GetLoc(uniform), static_cast<GLsizei>(mat.size()), false, (float*)&mat[0]); }
	inline void Uniform(const std::string& uniform, const std::vector<glm::mat4>& mat) { glUniformMatrix4fv(GetLoc(uniform), static_cast<GLsizei>(mat.size()), false, (float*)&mat[0]); }

	inline std::string GetName() { return name; }
	inline GLuint GetProgram() { return program; }
};

//
// MultipleFileShader
//
class MultipleFileShader : public Shader {
private:
	std::string vertexFileName;
	std::string fragmentFileName;
	std::string geometryFileName;
public:
	MultipleFileShader(const std::string& name);
	MultipleFileShader(const std::string& name, const std::string& vertexFileName, const std::string& fragmentFileName);
	MultipleFileShader(const std::string& name, const std::string& vertexFileName, const std::string& geometryFileName, const std::string& fragmentFileName);

	void Reload();
};

//
// ShaderManager
//
class ShaderManager {
private:
	std::unordered_map<std::string, Shader*> shaders;

	ShaderManager() {};
	ShaderManager(const ShaderManager&) = delete;
	ShaderManager& operator=(const ShaderManager&) = delete;
public:
	static ShaderManager& GetInstance() {
		static ShaderManager instance;
		return instance;
	}

	Shader* GetShader(const std::string& name);
	bool FindShader(const std::string& name);
	Shader* LoadShader(const std::string& name, const std::string& vertexFileName, const std::string& fragmentFileName);
	Shader* LoadShader(const std::string& name, const std::string& vertexFileName, const std::string& geometryFileName, const std::string& fragmentFileName);

	void Reload();

	Shader* CreateColorShader();
};



//
// ShaderBuffer
//
class ShaderBuffer {
private:
	GLuint buffer;
	GLenum target;

	void GenBuffer();
	void DeleteBuffer();

public:
	ShaderBuffer();
	~ShaderBuffer();

	void Bind(GLenum target);
	void Unbind(GLenum target);
	void BindBase(GLenum target, GLenum bindID);
	void WriteData(GLenum target, GLsizeiptr size, const GLvoid* data, GLenum usage);
	void ClearData(GLenum target, GLsizeiptr size, const GLvoid* data);
	void ReadData(GLenum target, GLsizeiptr size, GLvoid* data);
};
