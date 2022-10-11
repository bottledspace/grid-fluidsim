#define GLEW_STATIC
#include <GL/glew.h>
#include "Grid.hpp"
#include "GridRenderer.hpp"
#include "VectorFieldRenderer.hpp"
#include "SignedDistanceFieldRenderer.hpp"
#include <iostream>
#include <string>
#include <functional>
#include <chrono>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <GLFW/glfw3.h>

constexpr float kappa = 1e-6f;  // diffusion rate
constexpr float nu = 1e-6f; // kinematic viscosity
constexpr int N = 128-2;
constexpr float dt = 0.01f;
constexpr float h = 1.0f/N;
int width, height;
glm::ivec2 lastCursorPos;
glm::ivec2 cursorPos = {0,0};
Grid<glm::vec2> velocity(N,N);
Grid<float> density(N,N);
Grid<float> sdf(N,N,10000);

template <typename T>
void addSource(Grid<T>& dst, Grid<T>& src, float dt) {
    #pragma omp for collapse (2)
    for (int y = 0; y < N+2; y++)
    for (int x = 0; x < N+2; x++)
        dst(x,y) += src(x,y)*dt;
}

auto channel(Grid<glm::vec2>& grid, int dim) {
    Grid<float> res(grid.width(), grid.height());
    #pragma omp for collapse (2)
    for (int y = 0; y < grid.height()+2; y++)
    for (int x = 0; x < grid.width()+2; x++)
        res(x,y) = grid(x,y)[dim];
    return res;
}

void bounds(Grid<glm::vec2>& out) {
    const int w = out.width();
    const int h = out.height();
    #pragma omp for
    for (int y = 1; y <= h; y++) {
        out(  0,y) = glm::vec2(-1,1)*out(1,y);
        out(w+1,y) = glm::vec2(-1,1)*out(w,y);
    }
    #pragma omp for
    for (int x = 1; x <= w; x++) {
        out(x,  0) = glm::vec2(1,-1)*out(x,1);
        out(x,h+1) = glm::vec2(1,-1)*out(x,h);
    }
    out(0,0)     = 0.5f*(out(1,0)+out(0,1));
    out(w+1,0)   = 0.5f*(out(w,0)+out(w+1,1));
    out(w+1,h+1) = 0.5f*(out(w,h+1)+out(w+1,h));
    out(0,h+1)   = 0.5f*(out(1,h+1)+out(0,h));
}

void bounds(Grid<float>& out) {
    const int w = out.width();
    const int h = out.height();
    #pragma omp for
    for (int y = 1; y <= h; y++) {
        out(  0,y) = out(1,y);
        out(w+1,y) = out(w,y);
    }
    #pragma omp for
    for (int x = 1; x <= w; x++) {
        out(x,  0) = out(x,1);
        out(x,h+1) = out(x,h);
    }
    out(0,0)     = 0.5f*(out(1,0)+out(0,1));
    out(w+1,0)   = 0.5f*(out(w,0)+out(w+1,1));
    out(w+1,h+1) = 0.5f*(out(w,h+1)+out(w+1,h));
    out(0,h+1)   = 0.5f*(out(1,h+1)+out(0,h));
}

template <typename T>
auto diffuse(Grid<T>& grid, Grid<float>& sdf, float dt, float kappa) {
    const float a = dt*kappa*N*N;
    
    Grid<T> res(N,N);
    for (int k = 0; k < 20; k++) {
        #pragma for collapse (2)
        for (int y = 1; y <= N; y++)
        for (int x = 1; x <= N; x++) {
            if (sdf(x,y) > 0)
                continue;
            res(x,y) = (grid(x,y) +
                a*(res(x-1,y) + res(x+1,y)
                 + res(x,y-1) + res(x,y+1)))/(1+4*a);
        }
        bounds(res);
    }
    return res;
}

template <typename T>
auto advect(Grid<T>& grid, Grid<glm::vec2>& vel, float dt) {
    Grid<T> res(N,N);
    #pragma omp for collapse (2)
    for (int y = 1; y <= N; y++)
    for (int x = 1; x <= N; x++) {
        auto p = glm::vec2(x,y) - dt*N*vel(x,y);
        p.x = std::max(0.5f, std::min(N+0.5f, p.x));
        p.y = std::max(0.5f, std::min(N+0.5f, p.y));
        res(x,y) = grid(p);
    }
    bounds(res);
    return res;
}

// This version of advect restricts to the inner surface of a given sdf.
template <typename T>
auto advect(Grid<T>& grid, Grid<glm::vec2>& vel, Grid<float>& sdf, float dt) {
    Grid<T> res(N,N);
    #pragma omp for collapse (2)
    for (int y = 1; y <= N; y++)
    for (int x = 1; x <= N; x++) {
        auto p = glm::vec2(x,y) - dt*N*vel(x,y);
        p.x = std::max(0.5f, std::min(N+0.5f, p.x));
        p.y = std::max(0.5f, std::min(N+0.5f, p.y));
        if (sdf(p) <= 0)
            res(x,y) = grid(p);
    }
    bounds(res);
    return res;
}

auto project(Grid<glm::vec2>& grid, Grid<float>& sdf, float h) {    
    Grid<float> phi(grid.width(), grid.height());
    Grid<float> div(grid.width(), grid.height());
    Grid<glm::vec2> res(grid);
    
    // Solve ∇u = ΔΦ for Φ
    #pragma omp for collapse (2)
    for (int y = 1; y <= N; y++)
    for (int x = 1; x <= N; x++)
        if (sdf(x,y) <= 0)
        div(x,y) = -0.5f*h*(res(x+1,y).x-res(x-1,y).x
                           +res(x,y+1).y-res(x,y-1).y);
    bounds(div);
    for (int k = 0; k < 40; k++) {
        #pragma omp for collapse (2)
        for (int y = 1; y <= N; y++)
        for (int x = 1; x <= N; x++) {
            if (sdf(x,y) <= 0)
            phi(x,y) = (div(x,y)
                + phi(x+1,y) + phi(x-1,y)
                + phi(x,y+1) + phi(x,y-1))/4;
        }
        bounds(phi);
    }
    // Subtract ∇phi from res
    #pragma omp for collapse (2)
    for (int y = 1; y <= N; y++)
    for (int x = 1; x <= N; x++)
    if (sdf(x,y) <= 0) {
        res(x,y).x -= 0.5f*(phi(x+1,y)-phi(x-1,y))/h;
        res(x,y).y -= 0.5f*(phi(x,y+1)-phi(x,y-1))/h;
    }
    bounds(res);
    return res;
}

template <typename T>
auto extrapolate(Grid<T>& grid, Grid<float>& sdf) {
    Grid<int> valid0(N,N);
    Grid<int> valid1(N,N,0);
    Grid<T> res = grid;
    
    #pragma omp for collapse (2)
    for (int y = 0; y <= N+1; y++)
    for (int x = 0; x <= N+1; x++)
        valid0(x,y) = sdf(x,y) < 0;

    for (int k = 0; k < 100; k++) {
        #pragma omp for collapse (2)
        for (int y = 1; y <= N; y++)
        for (int x = 1; x <= N; x++) {
            if (valid0(x,y)) {
                valid1(x,y) = 1;
                continue;
            }
            T sum = T(0.0);
            float count = 0;
            if (x+1 <= N && valid0(x+1,y)) { sum += res(x+1,y); ++count; }
            if (x-1 >= 0 && valid0(x-1,y)) { sum += res(x-1,y); ++count; }
            if (y+1 <= N && valid0(x,y+1)) { sum += res(x,y+1); ++count; }
            if (y-1 >= 0 && valid0(x,y-1)) { sum += res(x,y-1); ++count; }
            if (count != 0) {
                res(x,y) = sum/count;
                valid1(x,y) = 1;
            }
        }
        swap(valid1, valid0);
    }
    return res;
}

auto adjustsdf(Grid<float>& sdf) {
    constexpr float T = 0.1*h;
    Grid<float> res = sdf;
    Grid<float> temp(sdf.width(), sdf.height());

    auto sq = [](float x) { return x*x; };
    
    // Integrate phi+sign(phi)(|grad(phi)|-1)=0
    for (int k = 0; k < 30; k++) {
        #pragma omp for collapse (2)
        for (int y = 1; y <= N; y++)
        for (int x = 1; x <= N; x++) {
            float phi = res(x,y);
            float dx1 = (phi - res(x-1,y))/h;
            float dx2 = (res(x+1,y) - phi)/h;
            float dy1 = (phi - res(x,y-1))/h;
            float dy2 = (res(x,y+1) - phi)/h;
            float s = phi/hypot(phi, h);
            
            phi -= T*std::max(s,0.0f)
                    *(sqrt(sq(std::max(dx1,0.0f))+sq(std::min(dx2,0.0f))
                          +sq(std::max(dy1,0.0f))+sq(std::min(dy2,0.0f))) - 1.0f);
            phi -= T*std::min(s,0.0f)
                    *(sqrt(sq(std::min(dx1,0.0f))+sq(std::max(dx2,0.0f))
                          +sq(std::min(dy1,0.0f))+sq(std::max(dy2,0.0f))) - 1.0f);
            temp(x,y) = phi;
        }
        bounds(res);
        swap(res, temp);
    }
    return res;
}

void simulate(Grid<glm::vec2>& velocity, Grid<float>& density, Grid<float>& sdf, float dt, float h) {
    
    // du/dt     = -(u.∇)u    + nu(∇^2u)  + g
    // variation = convection + diffusion + sources
    
    // Add gravity
    #pragma omp for collapse (2)
    for (int y = 1; y <= N; y++)
    for (int x = 1; x <= N; x++)
        if (sdf(x,y) <= 0)
        velocity(x,y).y -= 9.81*dt;
    bounds(velocity);
    
    //velocity = diffuse(velocity, sdf, dt, nu);
    //velocity = project(velocity, sdf, h);
    
    velocity = extrapolate(velocity, sdf);
    velocity = project(velocity, sdf, h);
    velocity = advect(velocity, velocity, sdf, dt);
    velocity = project(velocity, sdf, h);
    
    sdf = advect(sdf, velocity, dt);
    sdf = adjustsdf(sdf);
}

void fill_circle() {
    // Setup the SDF.
    auto dot2 = [](glm::vec2 a) { return dot(a,a); };
    auto sdHeart = [&](glm::vec2 p) {
        p.x = fabs(p.x);
        if( p.y+p.x>1.0f )
            return sqrtf(dot2(p-glm::vec2(0.25f,0.75f))) - sqrtf(2.0f)/4.0f;
        return sqrtf(fmin(dot2(p-glm::vec2(0.00f,1.00f)),
                          dot2(p-0.5f*fmax(p.x+p.y,0.0f)))) * copysign(1.0f,p.x-p.y);
    };
    
    constexpr float r = h*(N/4);
    for (int x = 0; x <= N+1; x++)
    for (int y = 0; y <= N+1; y++) {
        sdf(x,y) = fmin(sdf(x,y),sdHeart(2.0f*glm::vec2(x,y)/float(N)-glm::vec2(1.0f,0.5f)));
    }
}

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
    glfwGetFramebufferSize(window, &width, &height);
	glfwMakeContextCurrent(window);

	if (glewInit() != GLEW_OK) {
		std::cerr << "Failed to initialize GLEW." << std::endl;
		return EXIT_FAILURE;
	}

    GridRenderer densityRenderer;
    SignedDistanceFieldRenderer sdfRenderer;
    VectorFieldRenderer velocityRenderer;
    
    std::vector<uint8_t> pixels(width*height*3);
    FILE *pipeout = nullptr;
    if (argc == 2) {
        std::string command =
            "ffmpeg -y -f rawvideo -vcodec rawvideo -framerate 30"
            " -pix_fmt rgb24 -s 1024x1024 -i - -vf vflip -c:v libx264"
            " -crf 30 -b:v 0 -r 30 ";
        command += argv[1];  // Unsafe, should be escaped first...
        pipeout = popen(command.c_str(), "w");
        std::cout << width << "x" << height << std::endl;
    }

    fill_circle();

	glfwSwapInterval(1);
    int counter = 0;
	while (!glfwWindowShouldClose(window)) {
        if (++counter % 100 == 0) {
            fill_circle();
        }
        lastCursorPos = cursorPos;
    
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        cursorPos = glm::ivec2(N*xpos/512,N-N*ypos/512);
        cursorPos.x = std::max(0, std::min(N, cursorPos.x));
        cursorPos.y = std::max(0, std::min(N, cursorPos.y));
    
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT))
            velocity(cursorPos.x, cursorPos.y) += dt*5.0f*glm::vec2(cursorPos-lastCursorPos);
        
        simulate(velocity, density, sdf, dt, h);
        
		glClear(GL_COLOR_BUFFER_BIT);
        //densityRenderer.render(density);
        sdfRenderer.render(sdf);
        velocityRenderer.render(velocity);
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
