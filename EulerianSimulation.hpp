#pragma once

#include "Grid.hpp"
#include "GridRenderer.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>

class EulerianSimulation {
public:
    EulerianSimulation(int width, int height)
    : width(width),
      height(height),
      velocity(width, height),
      density(width, height)
    {}

    void update();
    void render() const;

    void fillCircle();
private:
    int width, height;
    Grid<glm::vec2> velocity;
    Grid<float> density;
    GridRenderer densityRenderer;
    float total;


    void diffuseVelocity();
    void advectVelocity();
    void advectDensity();
    void projectVelocity();
    void extrapolateVelocity();
    void addForces();
    void fixDensity();

    constexpr static float kappa = 1e-6f;  // diffusion rate
    constexpr static float nu = 1e-6f; // kinematic viscosity
    constexpr static float dt = 1.0f/30;
    constexpr static float h = 1.0f/128.0f;
};


float totalMass(Grid<float>& grid)
{
    float total = 0.0f;
    for (int y = 1; y <= grid.height(); y++)
    for (int x = 1; x <= grid.width(); x++)
        total += (grid(x,y)>0.5);
    return total;
}

void EulerianSimulation::fillCircle()
{
    float r = float(width)/6;
    for (int y = 0; y <= height+1; y++)
    for (int x = 0; x <= width+1; x++) {
        float dx = (width/2) - x;
        float dy = (height/3) - y;
        density(x,y) += (sqrt(dx*dx + dy*dy) < r);
    }
    total = totalMass(density);
}

void EulerianSimulation::diffuseVelocity()
{
    const float a = dt*kappa/(h*h);
    
    Grid<glm::vec2> res(width,height,glm::vec2(0.0f));
    for (int k = 0; k < 20; k++) {
        for (int y = 1; y <= height; y++)
        for (int x = 1; x <= width; x++) {
            res(x,y) = (velocity(x,y) +
                a*(res(x-1,y) + res(x+1,y)
                 + res(x,y-1) + res(x,y+1)))/(1+4*a);
        }
        bounds(res);
    }
    swap(velocity, res);
}

void EulerianSimulation::advectVelocity()
{
    Grid<glm::vec2> res(width,height);
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        auto p = glm::vec2(x,y) - dt*width*velocity(x,y);
        p.x = std::max(0.5f, std::min(width+0.5f, p.x));
        p.y = std::max(0.5f, std::min(height+0.5f, p.y));
        res(x,y) = velocity(p);
    }
    bounds(res);
    swap(res, velocity);
}

void EulerianSimulation::advectDensity()
{
    Grid<float> res(width,height,0);
    Grid<float> weights(width,height,0);

    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        auto p = glm::vec2(x,y) - dt*width*velocity(x,y);
        p.x = std::max(0.5f, std::min(width+0.5f, p.x));
        p.y = std::max(0.5f, std::min(height+0.5f, p.y));

        int x1 = (int)p.x;
        int y1 = (int)p.y;
        float u = p.x - x1;
        float v = p.y - y1;
        float w1 =        u*v;
        float w2 = (1.0f-u)*v;
        float w3 =        u*(1.0f-v);
        float w4 = (1.0f-u)*(1.0f-v);

        // Accumulate the amount we advect from each point
        weights(x1+1,y1+1) += w1;
        weights(x1,y1+1) += w2;
        weights(x1+1,y1) += w3;
        weights(x1,y1) += w4;
    }
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        auto p = glm::vec2(x,y) - dt*width*velocity(x,y);
        p.x = std::max(0.5f, std::min(width+0.5f, p.x));
        p.y = std::max(0.5f, std::min(height+0.5f, p.y));

        int x1 = (int)p.x;
        int y1 = (int)p.y;
        float u = p.x - x1;
        float v = p.y - y1;
        float w1 =        u*v;
        float w2 = (1.0f-u)*v;
        float w3 =        u*(1.0f-v);
        float w4 = (1.0f-u)*(1.0f-v);

        if (weights(x1+1,y1+1) > 1) w1 /= weights(x1+1,y1+1);
        if (weights(x1,y1+1) > 1)   w2 /= weights(x1,y1+1);
        if (weights(x1+1,y1) > 1)   w3 /= weights(x1+1,y1);
        if (weights(x1,y1) > 1)     w4 /= weights(x1,y1);
    
        res(x,y) += w1*density(x1+1,y1+1);
        res(x,y) += w2*density(x1,y1+1);
        res(x,y) += w3*density(x1+1,y1);
        res(x,y) += w4*density(x1,y1);
    }
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        float w = 1.0f-weights(x,y);
        if (w <= 0.0)
            continue;
        
        auto p = glm::vec2(x,y) + dt*width*velocity(x,y);
        p.x = std::max(0.5f, std::min(width+0.5f, p.x));
        p.y = std::max(0.5f, std::min(height+0.5f, p.y));

        int x1 = (int)p.x;
        int y1 = (int)p.y;
        float u = p.x - x1;
        float v = p.y - y1;
        float w1 =        u*v;
        float w2 = (1.0f-u)*v;
        float w3 =        u*(1.0f-v);
        float w4 = (1.0f-u)*(1.0f-v);

        res(x1+1,y1+1) += w*w1*density(x,y);
        res(x1,y1+1) += w*w2*density(x,y);
        res(x1+1,y1) += w*w3*density(x,y);
        res(x1,y1) += w*w4*density(x,y);
    }
    bounds(res);
    swap(res, density);
}

void EulerianSimulation::projectVelocity()
{    
    Grid<float> phi(width, height);
    Grid<float> div(width, height);
    Grid<glm::vec2> res = velocity;
    
    // Solve ∇u = ΔΦ for Φ
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++)
        if (density(x,y) >= 0.5f)
            div(x,y) = -0.5f*h*(res(x+1,y).x-res(x-1,y).x
                               +res(x,y+1).y-res(x,y-1).y);
    bounds(div);
    for (int k = 0; k < 50; k++) {
        for (int y = 1; y <= height; y++)
        for (int x = 1; x <= width; x++) {
            if (density(x,y) >= 0.5f)
                phi(x,y) = (div(x,y)
                    + phi(x+1,y) + phi(x-1,y)
                    + phi(x,y+1) + phi(x,y-1))/4;
        }
        bounds(phi);
    }
    // Subtract ∇phi from res
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++)
        if (density(x,y) >= 0.5f) {
            res(x,y).x -= 0.5f*(phi(x+1,y)-phi(x-1,y))/h;
            res(x,y).y -= 0.5f*(phi(x,y+1)-phi(x,y-1))/h;
        }
    bounds(res);
    swap(res, velocity);
}

void EulerianSimulation::extrapolateVelocity()
{
    Grid<int> valid0(width,height);
    Grid<int> valid1(width,height,0);
    Grid<glm::vec2> res = velocity;
    
    for (int y = 0; y <= height+1; y++)
    for (int x = 0; x <= width+1; x++)
        valid0(x,y) = density(x,y) >= 0.5;

    for (int k = 0; k < 50; k++) {
        for (int y = 0; y <= height+1; y++)
        for (int x = 0; x <= width+1; x++) {
            if (valid0(x,y)) {
                valid1(x,y) = 1;
                continue;
            }
            glm::vec2 sum(0.0f, 0.0f);
            float count = 0;
            if (x+1 < width && valid0(x+1,y))  sum += res(x+1,y), ++count;
            if (x-1 > 0 && valid0(x-1,y))      sum += res(x-1,y), ++count;
            if (y+1 < height && valid0(x,y+1)) sum += res(x,y+1), ++count;
            if (y-1 > 0 && valid0(x,y-1))      sum += res(x,y-1), ++count;
            if (count != 0) {
                res(x,y) = sum/count;
                valid1(x,y) = 1;
            }
        }
        swap(valid1, valid0);
    }
    bounds(res);
    swap(velocity, res);
}

void EulerianSimulation::addForces()
{
    // Add gravity
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++)
        if (density(x,y) >= 0.5)
            velocity(x,y).y -= 0.01;
    bounds(velocity);
}

void EulerianSimulation::fixDensity()
{
    float newTotalMass = totalMass(density);

    Grid<float> res(width,height);
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        float localSum = 0.0f;
        float totalWeight = 0.0f;
        for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            float weight = exp(-(dx*dx+dy*dy)/8.0f);
            localSum += weight*density(x+dx,y+dy);
            totalWeight += weight;
        }
        res(x,y) = localSum/totalWeight;
    }
    swap(res, density);
    for (int y = 1; y <= height; y++)
    for (int x = 1; x <= width; x++) {
        res(x,y) = 1.0/(1.0+exp(-(15.0*density(x,y)-7.5*newTotalMass/total)));
        /*float blur = (density(x,y) + density(x-1,y) + density(x+1,y)
            + density(x,y-1) + density(x,y+1))/5.0f;
        float alpha = 3.0f;
        res(x,y) = density(x,y) + alpha*(density(x,y) - blur);
        if (res(x,y) > 1)
            res(x,y) = 1;
        if (res(x,y) < 0)
            res(x,y) = 0;*/
    }
    bounds(res);
    swap(res, density);
}

void EulerianSimulation::update()
{
    addForces();
    advectVelocity();
    projectVelocity();
    extrapolateVelocity();

    advectDensity();
    fixDensity();
}

void EulerianSimulation::render() const
{
    densityRenderer.render(density);
}

#if 0
auto EulerianSimulation::adjustsdf(Grid<float>& sdf) {
    constexpr float T = 0.01*h;
    Grid<float> res = sdf;
    Grid<float> temp(sdf.width(), sdf.height());

    auto sq = [](float x) { return x*x; };
    
    // Integrate phi+sign(phi)(|grad(phi)|-1)=0
    for (int k = 0; k < 50; k++) {
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
#endif