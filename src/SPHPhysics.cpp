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
            const float *pos_y = particleSystem->getPositionsYData();
            for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
            {
                if (pos_y[i] > 490.0f)
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
                        const float *pos_x = particleSystem->getPositionsXData();
                        const float *pos_y = particleSystem->getPositionsYData();
                        std::vector<size_t> neighbors = grid->getNeighbors(pos_x[i], pos_y[i], 2 * h, particleSystem);

                        for (size_t neighborIdx : neighbors)
                        {
                            if (i == neighborIdx)
                                continue;
                            float dx = pos_x[i] - pos_x[neighborIdx];
                            float dy = pos_y[i] - pos_y[neighborIdx];
                            float distSqr = dx * dx + dy * dy;
                            if (distSqr < 4 * h2)
                                realNeighborCount++;
                        }

                        // Debug output commented out - vel data available via direct access if needed
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
        size_t particleCount = particleSystem->getParticleCount();

        // Use direct array pointers for maximum performance
        float *vel_x = particleSystem->getVelocitiesXData();
        float *vel_y = particleSystem->getVelocitiesYData();
        float *pos_x = particleSystem->getPositionsXData();
        float *pos_y = particleSystem->getPositionsYData();
        const float *accel_x = particleSystem->getAccelerationsXData();
        const float *accel_y = particleSystem->getAccelerationsYData();

        // Batch processing with SIMD optimization
        constexpr size_t BATCH_SIZE = 64;

#pragma omp parallel for schedule(static)
        for (size_t batch = 0; batch < particleCount; batch += BATCH_SIZE)
        {
            size_t batch_end = std::min(batch + BATCH_SIZE, particleCount);

// Vectorizable integration loop
#pragma omp simd aligned(vel_x, vel_y, pos_x, pos_y, accel_x, accel_y : 32)
            for (size_t i = batch; i < batch_end; ++i)
            {
                // Update velocity: v = v + a * dt
                vel_x[i] += accel_x[i] * dt;
                vel_y[i] += accel_y[i] * dt;

                // Update position: p = p + v * dt
                pos_x[i] += vel_x[i] * dt;
                pos_y[i] += vel_y[i] * dt;
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

        // Direct array access - eliminate expensive getters
        const std::vector<float> &pos_x = particleSystem->getPositionsX();
        const std::vector<float> &pos_y = particleSystem->getPositionsY();
        const std::vector<float> &masses = particleSystem->getMasses();

        float pos_i_x = pos_x[particleIndex];
        float pos_i_y = pos_y[particleIndex];
        float mass_i = masses[particleIndex];

        // Include self-contribution for density
        density += mass_i * POLY6_COEFF * h2 * h2 * h2;

        // Use cached neighbors for better performance
        const std::vector<size_t> &neighbors = particleSystem->getCachedNeighbors(particleIndex);

// Vectorizable loop with direct array access
#pragma omp simd reduction(+ : density)
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
            size_t neighborIdx = neighbors[j];
            if (particleIndex == neighborIdx)
                continue;

            float pos_j_x = pos_x[neighborIdx];
            float pos_j_y = pos_y[neighborIdx];
            float mass_j = masses[neighborIdx];

            float dx = pos_i_x - pos_j_x;
            float dy = pos_i_y - pos_j_y;
            float distSqr = dx * dx + dy * dy;

            if (distSqr < h2)
            {
                float h2_r2 = h2 - distSqr;
                density += mass_j * POLY6_COEFF * h2_r2 * h2_r2 * h2_r2;
            }
        }

        // Ensure minimum density
        density = std::max(density, 0.9f * restDensity);

        // Direct array access for setting values
        std::vector<float> &densities = particleSystem->getDensities();
        std::vector<float> &pressures = particleSystem->getPressures();

        densities[particleIndex] = density;

        float density_ratio = density / restDensity;
        // Basic EOS (Equation of State)
        float pressure_term = std::pow(density_ratio, gamma) - 1.0f;
        float pressure = (this->gasConstant * restDensity / gamma) * pressure_term;
        pressure = std::max(0.0f, pressure);
        pressures[particleIndex] = pressure;
    }

    void SPHPhysics::computeForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex)
    {
        // Kernel gradient constant (for Pressure AND Artificial Viscosity)
        const float SPIKY_GRAD_COEFF = -10.0f / (M_PI * std::pow(h, 5));
        const float EPS = 1e-6f;
        const float artificial_epsilon = 0.01f * h2;
        const float alpha = 1.0f;

        float pressureAccelX = 0.0f;
        float pressureAccelY = 0.0f;
        float viscosityAccelX = 0.0f;
        float viscosityAccelY = 0.0f;

        // Direct array access - eliminate expensive getters
        const std::vector<float> &pos_x = particleSystem->getPositionsX();
        const std::vector<float> &pos_y = particleSystem->getPositionsY();
        const std::vector<float> &vel_x = particleSystem->getVelocitiesX();
        const std::vector<float> &vel_y = particleSystem->getVelocitiesY();
        const std::vector<float> &densities = particleSystem->getDensities();
        const std::vector<float> &pressures = particleSystem->getPressures();
        const std::vector<float> &masses = particleSystem->getMasses();

        float pos_i_x = pos_x[particleIndex];
        float pos_i_y = pos_y[particleIndex];
        float vel_i_x = vel_x[particleIndex];
        float vel_i_y = vel_y[particleIndex];
        float density_i = densities[particleIndex];
        float pressure_i = pressures[particleIndex];

        // Make sure density isn't too close to zero before division
        if (density_i < EPS)
            density_i = EPS; // Safeguard

        // Use cached neighbors for better performance
        const std::vector<size_t> &neighbors = particleSystem->getCachedNeighbors(particleIndex);

        // Vectorizable loop with direct array access
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
            size_t neighborIdx = neighbors[j];
            if (particleIndex == neighborIdx)
                continue;

            // Direct array access for neighbor properties
            float pos_j_x = pos_x[neighborIdx];
            float pos_j_y = pos_y[neighborIdx];
            float vel_j_x = vel_x[neighborIdx];
            float vel_j_y = vel_y[neighborIdx];
            float density_j = densities[neighborIdx];
            float pressure_j = pressures[neighborIdx];
            float mass_j = masses[neighborIdx];

            float r_ij_x = pos_i_x - pos_j_x; // Vector from j to i
            float r_ij_y = pos_i_y - pos_j_y;
            float distSqr = r_ij_x * r_ij_x + r_ij_y * r_ij_y;

            if (distSqr < h2 && distSqr > EPS) // Check within h and avoid self/coincident
            {
                float dist = std::sqrt(distSqr);
                float dir_ij_x = r_ij_x / dist; // Normalized direction from j to i
                float dir_ij_y = r_ij_y / dist;

                // --- Pressure Acceleration Calculation ---
                float pressureTerm = -(pressure_i / (density_i * density_i) +
                                       pressure_j / (density_j * density_j));

                float h_r = h - dist;
                // Use Spiky gradient for pressure
                float gradW_spiky = SPIKY_GRAD_COEFF * h_r * h_r;
                pressureAccelX += mass_j * pressureTerm * gradW_spiky * dir_ij_x;
                pressureAccelY += mass_j * pressureTerm * gradW_spiky * dir_ij_y;

                // --- Monaghan Artificial Viscosity Acceleration ---
                float vel_ij_x = vel_i_x - vel_j_x; // Relative velocity v_i - v_j
                float vel_ij_y = vel_i_y - vel_j_y;
                float v_dot_r = vel_ij_x * r_ij_x + vel_ij_y * r_ij_y; // v_ij dot r_ij

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
                    viscosityAccelX += -mass_j * PI_ij * gradW_spiky * dir_ij_x;
                    viscosityAccelY += -mass_j * PI_ij * gradW_spiky * dir_ij_y;
                }
            }
        }

        // Set total acceleration directly in arrays
        std::vector<float> &accel_x = particleSystem->getAccelerationsX();
        std::vector<float> &accel_y = particleSystem->getAccelerationsY();

        accel_x[particleIndex] = pressureAccelX + viscosityAccelX + gravity.x;
        accel_y[particleIndex] = pressureAccelY + viscosityAccelY + gravity.y;
    }

    void SPHPhysics::computeBoundaryForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

        // Boundary parameters
        const float boundaryDistance = 1.5f * PARTICLE_RADIUS; // Detection distance
        const float boundaryStiffness = 10000.0f;
        const float boundaryDecay = 2.0f;

        // Direct array access
        const std::vector<float> &pos_x = particleSystem->getPositionsX();
        const std::vector<float> &pos_y = particleSystem->getPositionsY();

        float pos_x_val = pos_x[particleIndex];
        float pos_y_val = pos_y[particleIndex];

        float boundaryForceX = 0.0f;
        float boundaryForceY = 0.0f;

        // Left wall repulsion
        if (pos_x_val < boundaryDistance)
        {
            float dist = pos_x_val;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForceX += forceMagnitude;
        }

        // Right wall repulsion
        if (pos_x_val > width - boundaryDistance)
        {
            float dist = width - pos_x_val;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForceX -= forceMagnitude;
        }

        // Top wall repulsion
        if (pos_y_val < boundaryDistance)
        {
            float dist = pos_y_val;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForceY += forceMagnitude;
        }

        // Bottom wall repulsion
        if (pos_y_val > height - boundaryDistance)
        {
            float dist = height - pos_y_val;
            float normalizedDist = std::min(dist / boundaryDistance, 1.0f);
            float forceMagnitude = boundaryStiffness * std::pow(1.0f - normalizedDist, boundaryDecay);
            boundaryForceY -= forceMagnitude;
        }

        // Add boundary acceleration to existing acceleration using direct array access
        std::vector<float> &accel_x = particleSystem->getAccelerationsX();
        std::vector<float> &accel_y = particleSystem->getAccelerationsY();
        const std::vector<float> &densities = particleSystem->getDensities();

        float boundaryAccelX = boundaryForceX / densities[particleIndex];
        float boundaryAccelY = boundaryForceY / densities[particleIndex];

        accel_x[particleIndex] += boundaryAccelX;
        accel_y[particleIndex] += boundaryAccelY;
    }

    void SPHPhysics::integrateParticle(ParticleSystem *particleSystem, size_t particleIndex, float dt)
    {
        // Direct array access for integration - avoid getters/setters
        std::vector<float> &vel_x = particleSystem->getVelocitiesX();
        std::vector<float> &vel_y = particleSystem->getVelocitiesY();
        std::vector<float> &pos_x = particleSystem->getPositionsX();
        std::vector<float> &pos_y = particleSystem->getPositionsY();
        const std::vector<float> &accel_x = particleSystem->getAccelerationsX();
        const std::vector<float> &accel_y = particleSystem->getAccelerationsY();

        // Update velocity: v = v + a * dt
        vel_x[particleIndex] += accel_x[particleIndex] * dt;
        vel_y[particleIndex] += accel_y[particleIndex] * dt;

        // Update position: p = p + v * dt
        pos_x[particleIndex] += vel_x[particleIndex] * dt;
        pos_y[particleIndex] += vel_y[particleIndex] * dt;
    }

    void SPHPhysics::resolveCollisionsForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height)
    {
        constexpr float PARTICLE_RADIUS = 5.0f;

        // Direct array access
        std::vector<float> &pos_x = particleSystem->getPositionsX();
        std::vector<float> &pos_y = particleSystem->getPositionsY();
        std::vector<float> &vel_x = particleSystem->getVelocitiesX();
        std::vector<float> &vel_y = particleSystem->getVelocitiesY();

        float pos_x_val = pos_x[particleIndex];
        float pos_y_val = pos_y[particleIndex];
        float vel_x_val = vel_x[particleIndex];
        float vel_y_val = vel_y[particleIndex];

        // Simple boundary conditions with damping - failsafe only
        if (pos_x_val < PARTICLE_RADIUS)
        {
            pos_x_val = PARTICLE_RADIUS;
            vel_x_val = -vel_x_val * boundaryDamping;
        }
        if (pos_x_val > width - PARTICLE_RADIUS)
        {
            pos_x_val = width - PARTICLE_RADIUS;
            vel_x_val = -vel_x_val * boundaryDamping;
        }
        if (pos_y_val < PARTICLE_RADIUS)
        {
            pos_y_val = PARTICLE_RADIUS;
            vel_y_val = -vel_y_val * boundaryDamping;
        }
        if (pos_y_val > height - PARTICLE_RADIUS)
        {
            pos_y_val = height - PARTICLE_RADIUS;
            vel_y_val = -vel_y_val * boundaryDamping;
        }

        // Update arrays directly
        pos_x[particleIndex] = pos_x_val;
        pos_y[particleIndex] = pos_y_val;
        vel_x[particleIndex] = vel_x_val;
        vel_y[particleIndex] = vel_y_val;
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
