#define GLEW_STATIC
#include <GL/glew.h>
#include "EulerianSimulation.hpp"
#include <iostream>
#include <string>
#include <functional>
#include <chrono>
#include <cmath>
#include <GLFW/glfw3.h>


int main(int argc, char **argv)
{   
	if (!glfwInit()) {
		std::cerr << "Failed to initialize GLFW." << std::endl;
		return EXIT_FAILURE;
	}
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	GLFWwindow* window = glfwCreateWindow(512, 512, "grid", NULL, NULL);
	if (!window) {
		std::cerr << "Failed to create window." << std::endl;
		return EXIT_FAILURE;
	}
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
	glfwMakeContextCurrent(window);

	if (glewInit() != GLEW_OK) {
		std::cerr << "Failed to initialize GLEW." << std::endl;
		return EXIT_FAILURE;
	}
    
    std::vector<uint8_t> pixels(width*height*3);
    FILE *pipeout = nullptr;
    if (argc == 2) {
        std::string command =
            "ffmpeg -y -f rawvideo -vcodec rawvideo -framerate 30"
            " -pix_fmt rgb24 -s 512x512 -i - -vf vflip -c:v libx264"
            " -crf 30 -b:v 0 -r 30 ";
        command += argv[1];  // Unsafe, should be escaped first...
        pipeout = popen(command.c_str(), "w");
        std::cout << width << "x" << height << std::endl;
    }

    EulerianSimulation simulation(128-2,128-2);
    simulation.fillCircle();

	glfwSwapInterval(1);
    int counter = 0;
	while (!glfwWindowShouldClose(window)) {
        if (++counter == 500)
            simulation.fillCircle();
        simulation.update();
		glClear(GL_COLOR_BUFFER_BIT);
        simulation.render();
		glfwSwapBuffers(window);
        if (pipeout) {
            glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
            fwrite(pixels.data(), 1, pixels.size(), pipeout);
        }
		glfwPollEvents();
	}
    if (pipeout) {
        fflush(pipeout);
        pclose(pipeout);
    }
	glfwDestroyWindow(window);
	return 0;
}
