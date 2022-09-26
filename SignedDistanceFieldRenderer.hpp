#pragma once
#include "shader.hpp"

class SignedDistanceFieldRenderer {
public:
    SignedDistanceFieldRenderer()
    : m_prog(createShaderProgram(vertSrc, fragSrc)) {
    	// Create a Vertex Array Object to hold our mesh data.
    	glGenVertexArrays(1, &m_vao);
    	glBindVertexArray(m_vao);
    	// The VAO will point to a single Vertex Buffer Object which will interleave
    	// vertex information.
    	glGenBuffers(1, &m_vbo);
    	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    	// First: Spatial coordinates.
    	glEnableVertexAttribArray(0);
    	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0);
    	// Update the buffer
    	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);
        
        // Create a texture to send the field to the shader
    	glGenTextures(1, &m_texid);
        glBindTexture(GL_TEXTURE_2D, m_texid);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        // We hold off on completing the texture until the first render call
    }
    
    template <typename Grid>
    void render(const Grid& field) const {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_texid);
    	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F,
            field.width()+2, field.height()+2, 0,
    		GL_RED, GL_FLOAT, field.data());
        
		glUseProgram(m_prog);
        GLint loc = glGetUniformLocation(m_prog, "field");
		glUniform1i(loc, 0);
        glBindVertexArray(m_vao);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
private:
    GLuint m_prog;
    GLuint m_texid;
    GLuint m_vao;
    GLuint m_vbo;
    
    static const char* vertSrc;
    static const char* fragSrc;
    static const glm::vec2 quadVerts[4];
};
const char* SignedDistanceFieldRenderer::vertSrc = R"(#version 330
in vec2 uv;
out vec2 uvFrag;
void main() {
    uvFrag = uv;
	gl_Position = vec4(2.0*uv-1.0,0.0,1.0);
}
)";
const char* SignedDistanceFieldRenderer::fragSrc = R"(#version 330
uniform sampler2D field;
in vec2 uvFrag;
out vec3 color;
void main() {
    //float left = textureOffset(field, uvFrag, ivec2(-1,0)).r;
    //float right = textureOffset(field, uvFrag, ivec2(1,0)).r;
    //float up = textureOffset(field, uvFrag, ivec2(0,1)).r;
    //float down = textureOffset(field, uvFrag, ivec2(0,-1)).r;
	float value = texture(field, uvFrag).r;
    //if (left > 0 != right > 0 || up > 0 != down > 0)
    //   value += 0.5;
    color = (value > 0) ? 2.0*vec3(value, 0.0, 0.0) : 2.0*vec3(0.0, 0.0, -value);
}
)";
const glm::vec2 SignedDistanceFieldRenderer::quadVerts[4] = {
    {0.0f, 0.0f},
    {1.0f, 0.0f},
    {0.0f, 1.0f},
    {1.0f, 1.0f}
};