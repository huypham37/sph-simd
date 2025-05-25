#include "SPHPhysics.hpp"
#include "ParticleSystem.hpp"
#include "Grid.hpp"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <iostream>
#include <omp.h>

namespace sph
{
    SPHPhysics::SPHPhysics()
    {
        h = Config::SMOOTHING_RADIUS;
        h2 = h * h;
        viscosityCoefficient = Config::VISCOSITY;
        gasConstant = Config::GAS_CONSTANT;
        restDensity = Config::REST_DENSITY;
        boundaryDamping = Config::BOUNDARY_DAMPING;
        gamma = Config::GAMMA;
        timeStepCounter = 0;
    };

    void SPHPhysics::setSmoothingRadius(float smoothingRadius)
    {
        h = smoothingRadius;
        h2 = h * h;
    }

    void SPHPhysics::computeDensityPressure(ParticleSystem *particleSystem, Grid *grid)
    {
        // Increment time step counter at the beginning of each simulation step
        timeStepCounter++;
        bool shouldPrintDebug = (timeStepCounter % 1000 == 0);

        // Debug - select particles near the bottom of the screen
        debugParticles.clear();
        if (shouldPrintDebug && particleSystem->getParticleCount() > 0)
        {
            // Filter particles with y position greater than 490
            for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
            {
                Vector2f position = particleSystem->getPosition(i);

                if (position.y > 490.0f)
                {
                    debugParticles.insert(i);
                }
            }

            // Output how many particles we're debugging
            if (!debugParticles.empty())
            {
                std::cout << "Debugging " << debugParticles.size() << " particles near bottom" << std::endl;
            }
        }

// #pragma omp parallel
        {
            {
                int numThreads = omp_get_num_threads();
                if (timeStepCounter == 1 || timeStepCounter % 100 == 0)
                {
                    std::cout << "Running density/pressure calculation with " << numThreads << " threads" << std::endl;
                }
            }

            for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
            {
                computeDensityPressureForParticle(particleSystem, i);

                // Update the debug output code to respect the counter
                if (shouldPrintDebug && debugParticles.find(i) != debugParticles.end())
                {
                    {
                        int realNeighborCount = 0;
                        Vector2f pos = particleSystem->getPosition(i);
                        std::vector<size_t> neighbors = grid->getNeighbors(pos.x, pos.y, 2 * h, particleSystem);

                        for (size_t neighborIdx : neighbors)
                        {
                            if (i == neighborIdx)
                                continue;
                            Vector2f neighborPos = particleSystem->getPosition(neighborIdx);
                            Vector2f r = pos - neighborPos;
                            float distSqr = r.x * r.x + r.y * r.y;
                            if (distSqr < 4 * h2)
                                realNeighborCount++;
                        }

                        Vector2f vel = particleSystem->getVelocity(i);
                        // Debug output commented out
                    }
                }
            }
        }
    }

    void SPHPhysics::computeForces(ParticleSystem *particleSystem, Grid *grid)
    {
// #pragma omp simd
        for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
        {
            computeForcesForParticle(particleSystem, i);
        }
    }

    void SPHPhysics::integrate(ParticleSystem *particleSystem, float dt)
    {
        {
            {
                int numThreads = omp_get_num_threads();
                if (timeStepCounter == 1 || timeStepCounter % 500 == 0)
                {
                    std::cout << "Integrating with " << numThreads << " threads" << std::endl;
                }
            }

// #pragma omp simd
            for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
            {
                integrateParticle(particleSystem, i, dt);
            }
        }
    }

    void SPHPhysics::computeBoundaryForces(ParticleSystem *particleSystem, float width, float height)
    {
// #pragma omp simd
        for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
        {
            computeBoundaryForcesForParticle(particleSystem, i, width, height);
        }
    }

    void SPHPhysics::resolveCollisions(ParticleSystem *particleSystem, Grid *grid, float width, float height)
    {
// #pragma omp simd
        for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
        {
            resolveCollisionsForParticle(particleSystem, i, width, height);
        }
    }

    void SPHPhysics::computeDensityPressureForParticle(ParticleSystem *particleSystem, size_t particleIndex)
    {
        // Pre-compute kernel coefficient for optimization
        const float POLY6_COEFF = 4.0f / (M_PI * std::pow(h, 8));
        float density = 0.0f;

        // Get particle data
        Vector2f pos_i = particleSystem->getPosition(particleIndex);
        float mass_i = particleSystem->getMass(particleIndex);

        // Include self-contribution for density
        density += mass_i * POLY6_COEFF * h2 * h2 * h2;

        // Get neighbors and sum contribution
        std::vector<size_t> neighbors = particleSystem->getGrid()->getNeighbors(pos_i.x, pos_i.y, h, particleSystem);

// #pragma omp simd reduction(+ : density)
        for (size_t neighborIdx : neighbors)
        {
            if (particleIndex == neighborIdx)
                continue;

            Vector2f pos_j = particleSystem->getPosition(neighborIdx);
            float mass_j = particleSystem->getMass(neighborIdx);

            Vector2f r = pos_i - pos_j;
            float distSqr = r.x * r.x + r.y * r.y;

            if (distSqr < h2)
            {
                float h2_r2 = h2 - distSqr;
                density += mass_j * POLY6_COEFF * h2_r2 * h2_r2 * h2_r2;
            }
        }

        // Ensure minimum density
        density = std::max(density, 0.9f * restDensity);
        particleSystem->setDensity(particleIndex, density);

        float density_ratio = density / restDensity;
        // Basic EOS (Equation of State)
        float pressure_term = std::pow(density_ratio, gamma) - 1.0f;
        float pressure = (this->gasConstant * restDensity / gamma) * pressure_term;
        pressure = std::max(0.0f, pressure);
        particleSystem->setPressure(particleIndex, pressure);
    }

    void SPHPhysics::computeForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex)
    {
        // Kernel gradient constant (for Pressure AND Artificial Viscosity)
        const float SPIKY_GRAD_COEFF = -10.0f / (M_PI * std::pow(h, 5));
        const float EPS = 1e-6f;
        const float artificial_epsilon = 0.01f * h2;
        const float alpha = 1.0f;

        Vector2f pressureAcceleration(0.0f, 0.0f);
        Vector2f viscosityAcceleration(0.0f, 0.0f);

        // Get particle properties
        Vector2f pos_i = particleSystem->getPosition(particleIndex);
        Vector2f vel_i = particleSystem->getVelocity(particleIndex);
        float density_i = particleSystem->getDensity(particleIndex);
        float pressure_i = particleSystem->getPressure(particleIndex);

        // Make sure density isn't too close to zero before division
        if (density_i < EPS)
            density_i = EPS; // Safeguard

        // Get neighbors
        std::vector<size_t> neighbors = particleSystem->getGrid()->getNeighbors(pos_i.x, pos_i.y, h, particleSystem);

// #pragma omp simd
        for (size_t neighborIdx : neighbors)
        {
            if (particleIndex == neighborIdx)
                continue;

            // Get neighbor properties
            Vector2f pos_j = particleSystem->getPosition(neighborIdx);
            Vector2f vel_j = particleSystem->getVelocity(neighborIdx);
            float density_j = particleSystem->getDensity(neighborIdx);
            float pressure_j = particleSystem->getPressure(neighborIdx);
            float mass_j = particleSystem->getMass(neighborIdx);

            Vector2f r_ij = pos_i - pos_j; // Vector from j to i
            float distSqr = r_ij.x * r_ij.x + r_ij.y * r_ij.y;

            if (distSqr < h2 && distSqr > EPS) // Check within h and avoid self/coincident
            {
                float dist = std::sqrt(distSqr);
                Vector2f dir_ij = r_ij / dist; // Normalized direction from j to i

                // --- Pressure Acceleration Calculation ---
                float pressureTerm = -(pressure_i / (density_i * density_i) +
                                       pressure_j / (density_j * density_j));

                float h_r = h - dist;
                // Use Spiky gradient for pressure
                Vector2f gradW_spiky = SPIKY_GRAD_COEFF * h_r * h_r * dir_ij;
                pressureAcceleration += mass_j * pressureTerm * gradW_spiky;

                // --- Monaghan Artificial Viscosity Acceleration ---
                Vector2f vel_ij = vel_i - vel_j;                   // Relative velocity v_i - v_j
                float v_dot_r = vel_ij.x * r_ij.x + vel_ij.y * r_ij.y; // v_ij dot r_ij

                // Only apply viscosity if particles are approaching each other
                if (v_dot_r < 0.0f)
                {
                    // Average density
                    float rho_avg = 0.5f * (density_i + density_j);

                    // Speed of sound estimate
                    float c_i = std::sqrt(std::max(0.0f, pressure_i / density_i));
                    float c_j = std::sqrt(std::max(0.0f, pressure_j / density_j));
                    float c_avg = 0.5f * (c_i + c_j);

                    // Mu term from Monaghan's paper
                    float mu_ij = (h * v_dot_r) / (distSqr + artificial_epsilon);

                    // Viscosity term (PI_ij) - Ignoring beta term for liquids
                    float PI_ij = (-alpha * c_avg * mu_ij) / rho_avg;

                    // Viscosity Acceleration contribution
                    viscosityAcceleration += -mass_j * PI_ij * gradW_spiky;
                }
            }
        }

        // Set total acceleration = Pressure Acceleration + Viscosity Acceleration + Gravity
        particleSystem->setAcceleration(particleIndex, pressureAcceleration + viscosityAcceleration + gravity);
    }

    void SPHPhysics::computeBoundaryForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

        // Boundary parameters
        const float boundaryDistance = 1.5f * PARTICLE_RADIUS; // Detection distance
        const float boundaryStiffness = 10000.0f;
        const float boundaryDecay = 2.0f;

        Vector2f pos = particleSystem->getPosition(particleIndex);
        Vector2f boundaryForce(0.0f, 0.0f);

        // Left wall repulsion
        if (pos.x < boundaryDistance)
        {
            float dist = pos.x;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForce.x += forceMagnitude;
        }

        // Right wall repulsion
        if (pos.x > width - boundaryDistance)
        {
            float dist = width - pos.x;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForce.x -= forceMagnitude;
        }

        // Top wall repulsion
        if (pos.y < boundaryDistance)
        {
            float dist = pos.y;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForce.y += forceMagnitude;
        }

        // Bottom wall repulsion
        if (pos.y > height - boundaryDistance)
        {
            float dist = height - pos.y;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForce.y -= forceMagnitude;
        }

        // Add boundary acceleration to existing acceleration
        Vector2f currentAcc = particleSystem->getAcceleration(particleIndex);
        Vector2f boundaryAcc = boundaryForce / particleSystem->getDensity(particleIndex);
        particleSystem->setAcceleration(particleIndex, currentAcc + boundaryAcc);
    }

    void SPHPhysics::integrateParticle(ParticleSystem *particleSystem, size_t particleIndex, float dt)
    {
        Vector2f velocity = particleSystem->getVelocity(particleIndex) + particleSystem->getAcceleration(particleIndex) * dt;
        particleSystem->setVelocity(particleIndex, velocity);
        particleSystem->setPosition(particleIndex, particleSystem->getPosition(particleIndex) + velocity * dt);
    }

    void SPHPhysics::resolveCollisionsForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

        Vector2f pos = particleSystem->getPosition(particleIndex);
        Vector2f vel = particleSystem->getVelocity(particleIndex);

        // Simple boundary conditions with damping - failsafe only
        if (pos.x < PARTICLE_RADIUS)
        {
            pos.x = PARTICLE_RADIUS;
            vel.x = -vel.x * boundaryDamping;
        }
        if (pos.x > width - PARTICLE_RADIUS)
        {
            pos.x = width - PARTICLE_RADIUS;
            vel.x = -vel.x * boundaryDamping;
        }
        if (pos.y < PARTICLE_RADIUS)
        {
            pos.y = PARTICLE_RADIUS;
            vel.y = -vel.y * boundaryDamping;
        }
        if (pos.y > height - PARTICLE_RADIUS)
        {
            pos.y = height - PARTICLE_RADIUS;
            vel.y = -vel.y * boundaryDamping;
        }

        particleSystem->setPosition(particleIndex, pos);
        particleSystem->setVelocity(particleIndex, vel);
    }

    // Kernel functions
    float SPHPhysics::kernelPoly6(float distSquared)
    {
        if (distSquared >= h2)
            return 0.0f;

        const float POLY6_COEFF = 4.0f / (M_PI * std::pow(h, 8));
        float h2_r2 = h2 - distSquared;
        return POLY6_COEFF * h2_r2 * h2_r2 * h2_r2;
    }

    Vector2f SPHPhysics::kernelGradSpiky(float distance, const Vector2f &direction)
    {
        if (distance >= h || distance <= 0.0f)
            return Vector2f(0.0f, 0.0f);

        const float SPIKY_GRAD_COEFF = -10.0f / (M_PI * std::pow(h, 5));
        float h_r = h - distance;
        return SPIKY_GRAD_COEFF * h_r * h_r * direction;
    }

    float SPHPhysics::kernelViscosityLaplacian(float distance)
    {
        if (distance >= h)
            return 0.0f;

        const float VISC_LAP_COEFF = 40.0f / (M_PI * std::pow(h, 5));
        return VISC_LAP_COEFF * (h - distance);
    }
}
