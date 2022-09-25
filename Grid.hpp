#pragma once
#include <vector>
#include <glm/glm.hpp>

template <typename T>
class Grid {
public:
    Grid(const Grid &copy)
    : m_width(copy.m_width),
      m_height(copy.m_height),
      m_grid(copy.m_grid) {
    }
    Grid(Grid &&move)
    : m_width(std::exchange(move.m_width, 0)),
      m_height(std::exchange(move.m_height, 0)),
      m_grid(std::move(move.m_grid)) {
    }
    Grid(int width, int height, const T& val = T(0.0f))
    : m_width(width),
      m_height(height),
      m_grid((m_width+2)*(m_height+2), val) {
    }

    auto begin() { return m_grid.begin(); }
    auto end() { return m_grid.end(); }
    int width() const { return m_width; }
    int height() const { return m_height; }
    const T operator()(const glm::vec2& pt) const {
        int   x = (int)pt.x, y = (int)pt.y;
        float u = pt.x - x,  v = pt.y - y;
        return        u*v        * (*this)(x+1,y+1)
             + (1.0f-u)*v        * (*this)(x,y+1)
             +        u*(1.0f-v) * (*this)(x+1,y)
             + (1.0f-u)*(1.0f-v) * (*this)(x,y);
    }
    const T &operator()(int x, int y) const
        { return m_grid[y*(m_width+2)+x]; }
    T &operator()(int x, int y)
        { return m_grid[y*(m_width+2)+x]; }
    const T *data() const { return m_grid.data(); }
    T *data() { return m_grid.data(); }
    
    Grid &operator =(const Grid &copy) {
        m_width = copy.m_width;
        m_height = copy.m_height;
        m_grid = copy.m_grid;
        return *this;
    }
    
    friend void swap(Grid<T>& a, Grid<T>& b) {
        std::swap(a.m_width, b.m_width);
        std::swap(a.m_height, b.m_height);
        std::swap(a.m_grid, b.m_grid);
    }
private:
    int m_width, m_height;
    std::vector<T> m_grid;
};