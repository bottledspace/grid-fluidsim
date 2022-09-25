#pragma once
#include <iostream>

GLuint compileShader(GLenum type, const char* src) {
	// Compile shader
	GLuint sh = glCreateShader(type);
	glShaderSource(sh, 1, &src, NULL);
	glCompileShader(sh);

	// Handle errors and warnings
	GLint loglen, stat;
	glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &loglen);
	if (loglen > 0) {
		std::vector<char> buffer(loglen + 1);
		glGetShaderInfoLog(sh, loglen, NULL, buffer.data());
		std::cerr << buffer.data();
	}
	glGetShaderiv(sh, GL_COMPILE_STATUS, &stat);
	if (!stat) {
		std::cerr << "Failed to compile shader program." << std::endl;
		glDeleteShader(sh);
		return 0;
	}
	return sh;
}

GLuint createShaderProgram(const char* vertSrc, const char* fragSrc = nullptr) {
	GLuint prog, vert, frag;
	prog = glCreateProgram();
	vert = compileShader(GL_VERTEX_SHADER, vertSrc);
	if (fragSrc)
		frag = compileShader(GL_FRAGMENT_SHADER, fragSrc);
	glAttachShader(prog, vert);
	if (fragSrc)
		glAttachShader(prog, frag);
	glLinkProgram(prog);
	glDeleteShader(vert);
	if (fragSrc)
		glDeleteShader(frag);

	// Handle errors and warnings
	GLint loglen, stat;
	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &loglen);
	if (loglen > 0) {
		std::vector<char> buffer(loglen + 1);
		glGetProgramInfoLog(prog, loglen, NULL, buffer.data());
		std::cerr << buffer.data();
	}
	glGetProgramiv(prog, GL_LINK_STATUS, &stat);
	if (!stat) {
		std::cerr << "Failed to link shader program." << std::endl;
		glDeleteProgram(prog);
		return 0;
	}
	return prog;
}