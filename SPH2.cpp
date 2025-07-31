#include "SPH2.h"
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

const uint32_t TABLE_SIZE = 262144;


unsigned int myHash(vec2i cell)
{
    return ((uint)(cell.x * 73856093)
    ^ (uint)(cell.y * 19349663))
    % TABLE_SIZE;
}

float dist2(vec2 a, vec2 b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    return dx * dx + dy * dy;
}

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


bool operator!=(const vec2i& a, const vec2i& b)
{
    return a.x != b.x || a.y != b.y;
}










SPH::SPH()
{

    particleSquareWidth = 40;
    particleCount = particleSquareWidth * particleSquareWidth;
    //particles = new PARTICLE[particleCount];
    particleTable = new int[TABLE_SIZE];
    dirtTable = new int[TABLE_SIZE];


    // boxWidth = 0.5;
    // dt = 0.3;
    // mass = .02;
    // restDensity = 100;
    // gasConst = 1;
    // viscosity = 1044;
    // h = 0.001;
    // h2 = h * h;
    // g = -9.8;

    boxWidth = 200;
    dt = 0.01;
    mass = 1;
    restDensity = 1;
    gasConst = 20;
    viscosity = 0.5;
    h = 5;
    h2 = 25;
    g = -0.1;

    poly6 = 315.0 / (64.0 * M_PI * pow(h, 9));
    selfDens = mass * poly6 * pow(h, 6);
    spikyGrad = -45.0f / (M_PI * pow(h, 6));
    spikyLap = 45.0f / (M_PI * pow(h, 6));

    initParticles();
    createCircle(20);
}

SPH::~SPH()
{
    //delete[] particles;
    delete[] particleTable;
}


void SPH::add(float screen_x, float screen_y)
{
    static int count = 0;
    const int n = 1;
    const float dx = h;


    ++count;
    if (count < 8) return;
    count = 0;

    float x = (screen_x - 0.06) / 0.9 * boxWidth;
    float y = (screen_y - 0.02) / 0.9 * boxWidth;

    // x = 20;
    // y = 20;


    for (int i = -n; i <= n; ++i)
    {
        PARTICLE p;
        float offset = (((float) rand()) / float((RAND_MAX)) - 0.5) * dx * 0.1;
        p.position = {x + i * dx + offset, y - i * dx};
        p.velocity = {-15, -15};
        particles.push_back(p);
    }
}


void SPH::removeDirt(float screen_x, float screen_y) {
    float x = (screen_x - 0.06) / 0.9 * boxWidth;
    float y = (screen_y - 0.02) / 0.9 * boxWidth;

    PARTICLE p;
    p.position = { x, y };
    vec2i cell = cellOf(p);
    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            vec2i cellj = { cell.x + x, cell.y + y };
            
            unsigned cell_hash = myHash(cellj);

            
            int j = dirtTable[cell_hash];
            if (j == -1) continue;
            while (j < dirt.size())
            {
                
                PARTICLE& pj = dirt[j];
                if (pj.hash != cell_hash) break;
                float dist = dist2(p.position, pj.position);
                if (dist < h2)
                {
                    pj.removed = true;
                    removedDirt = true;
                }
                ++j;
            }

        }
    }
}

void SPH::update()
{
    cleanRemovedDirt();
    findNeighbors();
    calculateDensityAndPressure(particles);
    // calculateDensityAndPressure(dirt);
    calculateForces();
    updatePositions();
}

void SPH::render() {
    renderBalls();
    renderBeaker();
}

void SPH::renderBeaker() {

}

void SPH::renderBalls() {
    const float r = h / boxWidth * 0.9;
    //cout << "render" << endl;
    glBegin(GL_TRIANGLES);
//    glColor3f(116/256.0,204/256.0,244/256.0);
    glColor3f(35/256.0,137/256.0,218/256.0);
    for (int i = 0; i < particles.size(); ++i)
    {
        float x = particles[i].position.x;
        float y = particles[i].position.y;
        //cout << x << ", " << y << endl;
        float screenX = x / boxWidth * 0.9 + 0.06;
        float screenY = y / boxWidth * 0.9 + 0.02;
        for (int k = 0; k < num_triangles; ++k) {
            glVertex2f(screenX, screenY);
            glVertex2f(r * circle_coords[k].x + screenX, r * circle_coords[k].y + screenY);
            glVertex2f(r * circle_coords[k + 1].x + screenX, r * circle_coords[k + 1].y + screenY);
        }
    }
    glEnd();

    glBegin(GL_TRIANGLES);
    glColor3f(118/256.0,85/256.0,43/256.0);
    //glColor3f(1, 0, 0);
    for (int i = 0; i < dirt.size(); ++i)
    {
        float x = dirt[i].position.x;
        float y = dirt[i].position.y;
        
        float screenX = x / boxWidth * 0.9 + 0.06;
        float screenY = y / boxWidth * 0.9 + 0.02;
        for (int k = 0; k < num_triangles; ++k) {
            glVertex2f(screenX, screenY);
            glVertex2f(r * circle_coords[k].x + screenX, r * circle_coords[k].y + screenY);
            glVertex2f(r * circle_coords[k + 1].x + screenX, r * circle_coords[k + 1].y + screenY);
        }
    }
    glEnd();
}

//////////////////////////////////////////////////////
void SPH::cleanRemovedDirt() {
    if (!removedDirt) return;
    vector<PARTICLE> new_dirt;
    new_dirt.reserve(dirt.size());
    for (int i = 0; i < dirt.size(); ++i)
    {
        if (!dirt[i].removed)
        {
            new_dirt.push_back(dirt[i]);
        }
    }
    dirt = new_dirt;

    for (int i = 0; i < TABLE_SIZE; ++i)
    {
        dirtTable[i] = -1;
    }
    int last_hash = -1;
    for (int i = 0; i < dirt.size(); ++i)
    {
        unsigned int hash = dirt[i].hash;
        if (hash != last_hash)
        {
            dirtTable[hash] = i;
            last_hash = hash;
        }
    }
    removedDirt = false;
}

void SPH::createCircle(int n)
{
    float theta = 2 * M_PI / n;
    circle_coords = new vec2[n + 1];
    num_triangles = n;


    for (int i = 0; i < n + 1; ++i)
    {
        circle_coords[i] = { cos(i * theta), sin(i * theta) };
    }
}
void SPH::initParticles()
{
    particles.clear();
    dirt.clear();

    float dx = h / 3.0;
    int width = boxWidth / dx;
    int height = boxWidth / dx;

    float left = boxWidth * 2/ 7.0;
    float right = boxWidth * 5 / 7.0;

    float top = boxWidth * 4 / 5.0;
    float bottom = boxWidth / 5.0;

    float radius = boxWidth / 10.0;

    float wall_width = dx * 3;

    vec2 center = { boxWidth / 2, boxWidth / 2};
    vec2 bottomLeft = { left + radius, bottom + radius };
    vec2 bottomRight = { right - radius, bottom + radius };

    PARTICLE p;
    p.density = 0.0125335 * restDensity;
    p.pressure = gasConst * (p.density - restDensity);

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            // if (width / 5 < i && i < width * 4 / 5 &&
            //     height / 5 < j && j < height * 4 / 5)
            //     continue;
            p.position = {i * dx , (j + 5) * dx};
            float dist = sqrt(dist2(p.position, center));
            if (
                (abs(p.position.x - left) < wall_width || abs(p.position.x - right) < wall_width) && top - radius > p.position.y && p.position.y > bottom + radius ||
                (abs(p.position.y - bottom) < wall_width && left + radius < p.position.x && p.position.x < right - radius) ||
                (p.position.x < bottomLeft.x && p.position.y < bottomLeft.y && abs(sqrt(dist2(p.position, bottomLeft)) - radius) < wall_width) ||
                (p.position.x > bottomRight.x && p.position.y < bottomLeft.y && abs(sqrt(dist2(p.position, bottomRight)) - radius) < wall_width)
            ) {
            //if (dist < boxWidth * 0.4 || dist > boxWidth * 0.4 + dx * 3 ) continue;
                p.hash = myHash(cellOf(p));
                dirt.push_back(p);
            }
        }
    }

    // for (int i = 0; i < width; ++i) {
    //     for (int j = 0; j < height; ++j) {
    //         // if (width / 5 < i && i < width * 4 / 5 &&
    //         //     height / 5 < j && j < height * 4 / 5)
    //         //     continue;
    //         p.position = {i * dx , (j + 5) * dx};
    //         float dist = sqrt(dist2(p.position, center));
    //         if (dist < boxWidth * 0.4 || dist > boxWidth * 0.4 + dx * 3 ) continue;
    //         p.hash = myHash(cellOf(p));
    //         dirt.push_back(p);
    //     }
    // }

    sort(dirt.begin(), dirt.end(), [&](PARTICLE& l, PARTICLE&r)
    {
       return l.hash < r.hash; 
    });
    
    for (int i = 0; i < TABLE_SIZE; ++i)
    {
        dirtTable[i] = -1;
    }
    int last_hash = -1;
    for (int i = 0; i < dirt.size(); ++i)
    {
        unsigned int hash = dirt[i].hash;
        if (hash != last_hash)
        {
            dirtTable[hash] = i;
            last_hash = hash;
        }
    }

    
    //calculateDensityAndPressure(dirt);
}


vec2i SPH::cellOf(PARTICLE& p)
{
    return { (int)(p.position.x / h), (int)(p.position.y / h) };
}

void SPH::calculateHashes()
{
    float offset = h * 4;
    vector<PARTICLE> new_particles;
    for (int i = 0; i < particles.size(); ++i)
    {
        
        if (particles[i].position.x < -offset || particles[i].position.x > boxWidth + offset || particles[i].position.y < -offset)
            continue;
        PARTICLE& p = particles[i];
        vec2i cell = cellOf(p);
        p.hash = myHash(cell);
        new_particles.push_back(p);
    }

    particles = new_particles;
}

void SPH::sortParticles()
{
    sort(particles.begin(), particles.end(), [&](PARTICLE& l, PARTICLE&r)
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
        //dirtTable[i] = -1;
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

float SPH::W_poly6(float r2) {
    return poly6 * pow(h2 - r2, 3);
}


void SPH::calculateDensityAndPressure(std::vector<PARTICLE>& particles)
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
                        continue;
                    }
                    PARTICLE& pj = particles[j];
                    if (pj.hash != cell_hash) break;
                    float dist = dist2(pi.position, pj.position);
                    if (dist < h2)
                    {
                        density += mass * W_poly6(dist);
                    }
                    ++j;
                }

                /////////////////////
                j = dirtTable[cell_hash];
                if (j == -1) continue;
                while (j < dirt.size())
                {
                    PARTICLE& pj = dirt[j];
                    if (pj.hash != cell_hash) break;
                    float dist = dist2(pi.position, pj.position);
                    if (dist < h2)
                    {
                        density += mass * W_poly6(dist);
                    }
                    ++j;
                }
            }
        }
        pi.density = density + selfDens;
        // cout << pi.density / restDensity << endl;
        pi.pressure = gasConst * (pi.density - restDensity);
    }
}
void SPH::calculateForces()
{
    int count = 0;
    const float h2 = h * h;
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

                unsigned int hash = myHash(cellj);
                int j = particleTable[hash];
                if (j == -1) continue;

                while (j < particles.size())
                {
                    if (j == i)
                    {
                        ++j;
                        continue;
                    }
                    PARTICLE& pj = particles[j];
                    if (pj.hash != hash)
                    {
                        break;
                    }
                    float dist_2 = dist2(pi.position, pj.position);
                    if (dist_2 < h2)
                    {
                        float dist = sqrt(dist_2);
                        vec2 dir = (pj.position - pi.position).normalized();
                       

                        vec2 pressure_force = -dir * (mass * (pi.pressure + pj.pressure) / (2 * pj.density) * spikyGrad * pow(h - dist, 2));
                        pi.force = pi.force + pressure_force;

                        vec2 dv = pj.velocity - pi.velocity;
                        vec2 viscoForce = viscosity * mass * (dv / pj.density) * spikyLap * (h - dist);
                        pi.force = pi.force + viscoForce;
                        
                        ++count;
                        //vec2 surfaceForce = -sigma * neb_color * n / n.norm;
                    }
                    ++j;
                }
                ///////////////////////////////
                j = dirtTable[hash];
                if (j == -1) continue;

                while (j < dirt.size())
                {
                    PARTICLE& pj = dirt[j];
                    if (pj.hash != hash)
                    {
                        break;
                    }
                    float dist_2 = dist2(pi.position, pj.position);
                    if (dist_2 < h2)
                    {
                        float dist = sqrt(dist_2);
                        vec2 dir = (pj.position - pi.position).normalized();
                       

                        vec2 pressure_force = -dir * (mass * (pi.pressure + pj.pressure) / (2 * pj.density) * spikyGrad * pow(h - dist, 2));
                        pi.force = pi.force + pressure_force;

                        vec2 dv = pj.velocity - pi.velocity;
                        vec2 viscoForce = viscosity * mass * (dv / pj.density) * spikyLap * (h - dist);
                        pi.force = pi.force + viscoForce;

                       //cout << pressure_force.x << ", " << viscoForce.x << ", " << pressure_force.y << ", " << viscoForce.y << endl;
                        
                        ++count;
                        //vec2 surfaceForce = -sigma * neb_color * n / n.norm;
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
        //cout << p.density << endl;
        vec2 acceleration = (p.force + vec2(0, g)) / p.density;
        p.velocity = p.velocity + acceleration * dt;
        p.position = p.position + p.velocity * dt;

        // if (p.position.y < 0)
        // {
        //     p.position.y = 0;
        //     p.velocity.y = -p.velocity.y * elasticity;
        // }

        // if (p.position.y > boxWidth)
        // {
        //     p.position.y = boxWidth;
        //     p.velocity.y = -p.velocity.y * elasticity;
        // }

        // if (p.position.x < 0)
        // {
        //     p.position.x = 0;
        //     p.velocity.x = -p.velocity.x * elasticity;
        // }

        // if (p.position.x > boxWidth)
        // {
        //     p.position.x = boxWidth;
        //     p.velocity.x = -p.velocity.x * elasticity;
        // }
    }
}