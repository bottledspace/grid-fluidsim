#include "shader.hpp"

class VectorFieldRenderer {
public:
    VectorFieldRenderer()
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
    	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2)*2, NULL, GL_STATIC_DRAW);
        
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
    void render(const Grid& grid) const {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_texid);
    	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F,
            grid.width()+2, grid.height()+2, 0,
    		GL_RG, GL_FLOAT, grid.data());
        
		glUseProgram(m_prog);
        GLint loc = glGetUniformLocation(m_prog, "grid");
		glUniform1i(loc, 0);
        loc = glGetUniformLocation(m_prog, "width");
        glUniform1i(loc, 50);
        glBindVertexArray(m_vao);
		glDrawArraysInstanced(GL_LINES, 0, 2, 50*50);
    }
private:
    GLuint m_prog;
    GLuint m_texid;
    GLuint m_vao;
    GLuint m_vbo;
    
    static const char* vertSrc;
    static const char* fragSrc;
};
const char* VectorFieldRenderer::vertSrc = R"(#version 330
uniform int width;
uniform sampler2D grid;
void main() {
    float xi = float(gl_InstanceID % width);
    float yi = float(gl_InstanceID / width);
    vec2 uv = vec2(xi,yi)/width;
    vec2 dir = 2.0*texture(grid, uv).xy/width;
    vec2 pos = 2.0*uv-1.0;
	gl_Position = vec4(pos + gl_VertexID*dir,0.0,1.0);
}
)";
const char* VectorFieldRenderer::fragSrc = R"(#version 330
out vec3 color;
void main() {
    color = vec3(0.0,0.0,0.0);
}
)";