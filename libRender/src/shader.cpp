#include "shader.h"

//
// Shader
//
Shader::Shader(const std::string& name) : name(name) {}

void Shader::Use() {
	glUseProgram(program);
}
void Shader::Unuse() {
	glUseProgram(0);
}

void Shader::CompileShader(const std::string& vertexShaderCode, const std::string& fragmentShaderCode) {
	GLuint vertexShaderID = CompileAndCheck(vertexShaderCode, GL_VERTEX_SHADER);
	GLuint fragmentShaderID = CompileAndCheck(fragmentShaderCode, GL_FRAGMENT_SHADER);
	GLuint programID = LinkAndCheck(vertexShaderID, fragmentShaderID);
	program = programID;
}

void Shader::CompileShader(const std::string& vertexShaderCode, const std::string& geometryShaderCode, const std::string& fragmentShaderCode) {
	GLuint vertexShaderID = CompileAndCheck(vertexShaderCode, GL_VERTEX_SHADER);
	GLuint geometryShaderID = CompileAndCheck(geometryShaderCode, GL_GEOMETRY_SHADER);
	GLuint fragmentShaderID = CompileAndCheck(fragmentShaderCode, GL_FRAGMENT_SHADER);
	GLuint programID = LinkAndCheck(vertexShaderID, geometryShaderID, fragmentShaderID);
	program = programID;
}

GLuint Shader::CompileAndCheck(const std::string& shaderCode, GLenum shaderType) {
	GLuint shaderID = glCreateShader(shaderType);

	GLint result = GL_FALSE;
	int infoLogLength;

	// Compile Vertex Shader
	char const* SourcePointer = shaderCode.c_str();
	glShaderSource(shaderID, 1, &SourcePointer, NULL);
	glCompileShader(shaderID);

	// Check Vertex Shader
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if (infoLogLength > 1) {
		GLchar* errorMsg = new GLchar[infoLogLength];
		glGetShaderInfoLog(shaderID, infoLogLength, NULL, &errorMsg[0]);
		PRINT_ERR(errorMsg);
	} else {
		//PRINT("shader compiled");
	}
	return shaderID;
}

GLuint Shader::LinkAndCheck(GLuint vertexShaderID, GLuint fragmentShaderID) {
	GLint result = GL_FALSE;
	int infoLogLength;

	// Link the program
	GLuint programID = glCreateProgram();
	glAttachShader(programID, vertexShaderID);
	glAttachShader(programID, fragmentShaderID);
	glLinkProgram(programID);

	// Check the program
	glGetProgramiv(programID, GL_LINK_STATUS, &result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if (infoLogLength > 1) {
		GLchar* errorMsg = new GLchar[infoLogLength];
		glGetProgramInfoLog(programID, infoLogLength, NULL, &errorMsg[0]);
		PRINT_ERR(errorMsg);
	} else {
		//PRINT("program linked");
	}

	glDeleteShader(vertexShaderID);
	glDeleteShader(fragmentShaderID);

	return programID;
}

GLuint Shader::LinkAndCheck(GLuint vertexShaderID, GLuint geometryShaderID, GLuint fragmentShaderID) {
	GLint result = GL_FALSE;
	int infoLogLength;

	// Link the program
	GLuint programID = glCreateProgram();
	glAttachShader(programID, vertexShaderID);
	glAttachShader(programID, geometryShaderID);
	glAttachShader(programID, fragmentShaderID);
	glLinkProgram(programID);

	// Check the program
	glGetProgramiv(programID, GL_LINK_STATUS, &result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if (infoLogLength > 1) {
		GLchar* errorMsg = new GLchar[infoLogLength];
		glGetProgramInfoLog(programID, infoLogLength, NULL, &errorMsg[0]);
		PRINT_ERR(errorMsg);
	} else {
		//PRINT("program linked");
	}

	glDeleteShader(vertexShaderID);
	glDeleteShader(geometryShaderID);
	glDeleteShader(fragmentShaderID);

	return programID;
}

//
// MultipleFileShader
//
MultipleFileShader::MultipleFileShader(const std::string& name) :
	Shader(name), vertexFileName(""), fragmentFileName(""), geometryFileName("") {}

MultipleFileShader::MultipleFileShader(const std::string& name, const std::string& vertexFileName, const std::string& fragmentFileName) :
	Shader(name), vertexFileName(vertexFileName), fragmentFileName(fragmentFileName), geometryFileName("") {
	Reload();
}

MultipleFileShader::MultipleFileShader(const std::string& name, const std::string& vertexFileName, const std::string& geometryFileName, const std::string& fragmentFileName) :
	Shader(name), vertexFileName(vertexFileName), fragmentFileName(fragmentFileName), geometryFileName(geometryFileName) {
	Reload();
}

void MultipleFileShader::Reload() {
	if (vertexFileName == "" || fragmentFileName == "") { // do not reload generated shaders
		return;
	}

	PRINT("Loading shader: " + name);

	std::string vertexShaderCode = ReadFromFile(vertexFileName);
	std::string fragmentShaderCode = ReadFromFile(fragmentFileName);

	if (geometryFileName != "") {
		std::string geometryShaderCode = ReadFromFile(geometryFileName);
		CompileShader(vertexShaderCode, geometryShaderCode, fragmentShaderCode);
	} else {
		CompileShader(vertexShaderCode, fragmentShaderCode);
	}

}


//
// ShaderManager
//
Shader* ShaderManager::GetShader(const std::string& name) {
	if (!FindShader(name)) {
		PRINT_ERR("Shader does not exist: " + name);
		return nullptr;
	} else {
		return shaders.at(name);
	}
}

bool ShaderManager::FindShader(const std::string& name) {
	if (shaders.find(name) == shaders.end()) {
		return false;
	} else {
		return true;
	}
}

Shader* ShaderManager::LoadShader(const std::string& name, const std::string& vertexFileName, const std::string& fragmentFileName) {
	return LoadShader(name, vertexFileName, "", fragmentFileName);
}
Shader* ShaderManager::LoadShader(const std::string& name, const std::string& vertexFileName, const std::string& geometryFileName, const std::string& fragmentFileName) {
	if (!FindShader(name)) {
		shaders[name] = new MultipleFileShader(name, vertexFileName, geometryFileName, fragmentFileName);
		return shaders[name];
	} else {
		PRINT_ERR("Shader already exists: " + name);
		return nullptr;
	}
}

void ShaderManager::Reload() {
	for (auto& entry : shaders) {
		entry.second->Reload();
	}
}


Shader* ShaderManager::CreateColorShader() {
	const std::string name = "genericColorShader";

	if (FindShader(name)) {
		return shaders.at(name);
	}

	PRINT("Creating shader: " + name);
	Shader* s = new MultipleFileShader(name);

	std::string vertexShaderCode = "#version 330 core\n"
		"layout(location = 0) in vec3 in_position;\n"
		"uniform mat4 MVP;\n"
		"void main() {\n"
		"gl_Position = MVP * vec4(in_position, 1.0);\n"
		"}";

	std::string fragmentShaderCode = "#version 330 core\n"
		"out vec4 out_col;\n"
		"uniform vec3 color;\n"
		"void main(void) {\n"
		"out_col = vec4(color, 1.0);\n"
		"}";

	s->CompileShader(vertexShaderCode, fragmentShaderCode);
	shaders[name] = s;

	return s;
}


//
// ShaderBuffer
//
ShaderBuffer::ShaderBuffer() {
	GenBuffer();
}
ShaderBuffer::~ShaderBuffer() {
	DeleteBuffer();
}
void ShaderBuffer::GenBuffer() {
	glGenBuffers(1, &buffer);
}
void ShaderBuffer::DeleteBuffer() {
	glDeleteBuffers(1, &buffer);
}
void ShaderBuffer::Bind(GLenum target) {
	glBindBuffer(target, buffer);
}
void ShaderBuffer::Unbind(GLenum target) {
	glBindBuffer(target, 0);
}
void ShaderBuffer::BindBase(GLenum target, GLenum bindID) {
	glBindBufferBase(target, bindID, buffer);
}
// target = GL_UNIFORM_BUFFER | GL_SHADER_STORAGE_BUFFER | GL_DRAW_INDIRECT_BUFFER
// usage = GL_DYNAMIC_DRAW | GL_STATIC_DRAW
void ShaderBuffer::WriteData(GLenum target, GLsizeiptr size, const GLvoid* data, GLenum usage) {
	Bind(target);
	glBufferData(target, size, data, usage);
	Unbind(target);
}
void ShaderBuffer::ClearData(GLenum target, GLsizeiptr size, const GLvoid* data) {
	Bind(target);
	glBufferSubData(target, 0, size, data);
	Unbind(target);
}
void ShaderBuffer::ReadData(GLenum target, GLsizeiptr size, GLvoid* data) {
	Bind(target);
	glGetBufferSubData(target, 0, size, data);
	Unbind(target);
}
