#include "SPH.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

vec2 operator+(const vec2& a, const vec2& b)
{
    return { a.x + b.x, a.y + b.y };
}

vec2 operator-(const vec2& a, const vec2& b)
{
    return { a.x - b.x, a.y - b.y };
}

vec2 operator-(const vec2& a)
{
    return { -a.x, -a.y };
}

vec2 operator*(const vec2& a, float k)
{
    return { k * a.x, k * a.y };
}

vec2 operator*(float k, const vec2& a)
{
    return { k * a.x, k * a.y };
}

vec2 operator/(const vec2& a, float k)
{
    return { a.x/k, a.y/k };
}

vec2i operator+(const vec2i& a, const vec2i& b)
{
    return { a.x + b.x, a.y + b.y };
}












const uint32_t TABLE_SIZE = 262144;
SPH::SPH()
{
    boxWidth = 0.5;

    particleSquareWidth = 100;
    particleCount = particleSquareWidth * particleSquareWidth;
    //particles = new PARTICLE[particleCount];
    particleTable = new int[TABLE_SIZE];

    //sphere = new Geometry("resources/lowsphere.obj");
    //sphereModelMtxs = new glm::mat4[particleCount];


    dt = 0.003;
    mass = .02;
    restDensity = 1000;
    gasConst = 0.1;
    viscosity = 1.044;
    h = boxWidth / 100.0;
    g = -9.8;
    tension = 0.2;

    
    // boxWidth = 100;
    // dt = 0.01;
    // mass = 1;
    // restDensity = 1;
    // gasConst = 20;
    // viscosity = 0.5;
    // h = 5;
    // g = -0.1;

    poly6 = 315.0 / (64.0 * M_PI * pow(h, 9));
    selfDens = mass * poly6 * pow(h, 6);
    spikyGrad = -45.0f / (M_PI * pow(h, 6));
    spikyLap = 45.0f / (M_PI * pow(h, 6));

    initParticles();
}

SPH::~SPH()
{
    //delete[] particles;
    delete[] particleTable;
}

void SPH::add(float screen_x, float screen_y)
{
    static int count = 0;
    const int n = 3;
    const float dx = h;

    ++count;
    if (count < 4) return;
    count = 0;

    float x = (screen_x - 0.06) / 0.9 * boxWidth;
    float y = (screen_y - 0.06) / 0.9 * boxWidth;


    for (int i = -n; i <= n; ++i)
    {
        PARTICLE p;
        float offset = (((float) rand()) / float((RAND_MAX)) - 0.5) * dx * 0.1;
        p.position = {x + i * dx + offset, y};
        p.velocity = {0, -0.4};
        particles.push_back(p);
    }
}

void SPH::update()
{
    findNeighbors();
    calculateDensityAndPressure();
    calculateForces();
    updatePositions();
}

void SPH::render() {

    glBegin(GL_QUADS);
    glColor3f(1, 1, 1);
    for (int i = 0; i < particles.size(); ++i)
    {
        const float r = 0.001;
        float x = particles[i].position.x;
        float y = particles[i].position.y;
        //cout << x << ", " << y << endl;
        float screenX = x / boxWidth * 0.9 + 0.06;
        float screenY = y / boxWidth * 0.9 + 0.06;
        glVertex2f(screenX - r, screenY - r);
        glVertex2f(screenX + r, screenY - r);
        glVertex2f(screenX + r, screenY + r);
        glVertex2f(screenX - r, screenY + r);
    }
    glEnd();
}

//////////////////////////////////////////////////////

void SPH::initParticles()
{
    particles.clear();
    srand(1000);
    float particleSeparation = boxWidth / particleSquareWidth;
    for (int i = 0; i < particleSquareWidth; ++i)
    {
        for (int j = 0; j < particleSquareWidth; ++j)
        {
            float x = i * particleSeparation + (((float) rand()) / float((RAND_MAX)) * 0.5 - 1) * particleSeparation / 5.0;
            float y = j * particleSeparation + (((float) rand()) / float((RAND_MAX)) * 0.5 - 1) * particleSeparation / 5.0;
            
            PARTICLE p;
            p.position = {x / 4 + 0.375, y / 4 + 0.375};
            p.velocity = {0, 0};

            //particles.push_back(p);
        }
    }
}

unsigned int myHash(vec2i cell)
{
    return ((uint)(cell.x * 73856093)
    ^ (uint)(cell.y * 19349663))
    % TABLE_SIZE;
}

vec2i SPH::cellOf(PARTICLE& p)
{
    return { (int)(p.position.x / h), (int)(p.position.y / h) };
}

void SPH::calculateHashes()
{
    for (int i = 0; i < particles.size(); ++i)
    {
        PARTICLE& p = particles[i];
        vec2i cell = cellOf(p);
        p.hash = myHash(cell);
    }
}

void SPH::sortParticles()
{
    sort(particles.begin(), particles.end(), [](PARTICLE& l, PARTICLE&r)
    {
       return l.hash < r.hash; 
    });
}

void SPH::findNeighbors()
{
    calculateHashes();
    sortParticles();
    for (int i = 0; i < TABLE_SIZE; ++i)
    {
        particleTable[i] = -1;
    }

    unsigned int last_hash = -1;
    for (int i = 0; i < particles.size(); ++i)
    {
        unsigned int hash = particles[i].hash;
        if (hash != last_hash)
        {
            particleTable[hash] = i;
            last_hash = hash;
        }
    }
}

float dist2(vec2 a, vec2 b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    return dx * dx + dy * dy;
}

void SPH::calculateDensityAndPressure()
{
    float h2 = h * h;
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i)
    {
        float density = 0;
        PARTICLE& pi = particles[i];
        vec2i cell = cellOf(pi);
        for (int x = -1; x <= 1; ++x)
        {
            for (int y = -1; y <= 1; ++y)
            {
                vec2i cellj = { cell.x + x, cell.y + y };
                
                unsigned cell_hash = myHash(cellj);

                
                int j = particleTable[cell_hash];
                if (j == -1) continue;
                while (j < particles.size())
                {
                    if (j == i)
                    {
                        ++j;
                        //continue;
                    }
                    PARTICLE& pj = particles[j];
                    if (pj.hash != cell_hash) break;
                    float dist = dist2(pi.position, pj.position);
                    if (dist < h2)
                    {
                        density += mass * poly6 * pow(h2 - dist, 3);
                    }
                    ++j;
                }
            }
        }
        pi.density = density + selfDens;
        pi.pressure = gasConst * (pi.density - restDensity);
    }
}

void SPH::calculateForces()
{
    float h2 = h * h;
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i)
    {
        PARTICLE& pi = particles[i];
        pi.force = { 0, 0 };
        vec2i cell = cellOf(pi);

        for (int x = -1; x <= 1; ++x)
        {
            for (int y = -1; y <= 1; ++y)
            {
                vec2i cellj = { cell.x + x, cell.y + y };
                // if (cellj.y == -1)
                // {
                //     float dist = h / 2;
                //     vec2 dir = { 1.0, 0.0 };

                //     vec2 pressure_force = -dir * (mass * (pi.pressure + 0) / (2 * restDensity) * spikyGrad)* pow(h - dist, 2);
                //     pi.force = pi.force + pressure_force;

                //     vec2 dv = - pi.velocity;
                //     vec2 viscoForce = viscosity * mass * (dv / restDensity) * spikyLap * (h - dist);
                //     pi.force = pi.force + viscoForce;
                // }
                unsigned int hash = myHash(cellj);
                int j = particleTable[hash];
                if (j == -1) continue;

                while (j < particles.size())
                {
                    if (j == i)
                    {
                        ++j;
                        //continue;
                    }
                    PARTICLE& pj = particles[j];
                    if (pj.hash != pi.hash)
                    {
                        break;
                    }
                    float dist = dist2(pi.position, pj.position);
                    if (dist < h2)
                    {
                        dist = sqrt(dist);
                        vec2 dir = (pi.position - pj.position).normalized();

                        vec2 pressure_force = -dir * (mass * (pi.pressure + pj.pressure) / (2 * pj.density) * spikyGrad)* pow(h - dist, 2);
                        pi.force = pi.force + pressure_force;

                        vec2 dv = pj.velocity - pi.velocity;
                        vec2 viscoForce = viscosity * mass * (dv / pj.density) * spikyLap * (h - dist);
                        pi.force = pi.force + viscoForce;
                        
                    }
                    ++j;
                }
            }
        }
    }
}

void SPH::updatePositions()
{
    float elasticity = 0.5f;
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i)
    {
        PARTICLE& p = particles[i];
        //cout << p.force.x << ", " << p.force.y << endl;
        cout << p.density << endl;
        vec2 acceleration = (p.force) / p.density + vec2(0, g);
        p.velocity = p.velocity + acceleration * dt;
        p.position = p.position + p.velocity * dt;

        if (p.position.y < 0)
        {
            p.position.y = 0;
            p.velocity.y = -p.velocity.y * elasticity;
        }

        if (p.position.y > boxWidth)
        {
            p.position.y = boxWidth;
            p.velocity.y = -p.velocity.y * elasticity;
        }

        if (p.position.x < 0)
        {
            p.position.x = 0;
            p.velocity.x = -p.velocity.x * elasticity;
        }

        if (p.position.x > boxWidth)
        {
            p.position.x = boxWidth;
            p.velocity.x = -p.velocity.x * elasticity;
        }
    }
}