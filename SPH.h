#pragma once
#include <cmath>
#include <vector>


struct vec2
{
    float x;
    float y;

    vec2() : x(0), y(0) {}
    vec2(float x, float y) : x(x), y(y) {}
    float norm()
    {
        return sqrt(x * x + y * y);
    }

    vec2 normalized()
    {
        float norm = this->norm();
        if (norm == 0) return { 0, 0 };
        else return { x / norm, y / norm };
    }
};

struct vec2i
{
    int x;
    int y;
};

vec2 operator+(const vec2& a, const vec2& b);
vec2 operator-(const vec2& a, const vec2& b);
vec2 operator-(const vec2& a);
vec2 operator*(const vec2& a, float k);
vec2 operator*(float k, const vec2& a);
vec2 operator/(const vec2& a, float k);
vec2i operator+(const vec2i& a, const vec2i& b);


struct PARTICLE
{
    vec2 position;
    vec2 velocity;
    vec2 acceleration;
    vec2 force;

    float density;
    float pressure;
    unsigned int hash;
};

class SPH
{
public:
    SPH();
    ~SPH();
    void add(float x, float y);
    void reset() { initParticles(); }
    void update();
    void render();
private:
    vec2i cellOf(PARTICLE&);

    void initParticles();
    void calculateHashes();
    void sortParticles();
    void findNeighbors();
    void calculateDensityAndPressure();
    void calculateForces();
    void updatePositions();


    float boxWidth;

    int particleSquareWidth;
    int particleCount;
    std::vector<PARTICLE> particles;
    int* particleTable;

    float dt;
    float mass;
    float restDensity;
    float gasConst;
    float viscosity;
    float h;
    float g;
    float tension;

    float poly6;
    float selfDens;
    float spikyGrad;
    float spikyLap;
};