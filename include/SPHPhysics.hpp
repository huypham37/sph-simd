#pragma once

#include "Grid.hpp" // Added missing include for Grid class and Vector2f
#include "SPHConfig.hpp"
#include <unordered_set> // Added missing include for unordered_set
#include <vector>

namespace sph
{

  // Forward declarations
  class ParticleSystem;

  /**
   * @brief Handles SPH physics calculations using Structure of Arrays
   *
   * Responsible for implementing the SPH fluid simulation algorithm,
   * including density, pressure, forces, and physical interactions
   * optimized for SIMD operations with SoA layout.
   */
  class SPHPhysics
  {
  public:
    /**
     * @brief Constructor
     */
    SPHPhysics();
    /**
     * @brief Set simulation parameters
     *
     * @param smoothingRadius SPH kernel smoothing radius
     */
    void setSmoothingRadius(float smoothingRadius);

    /**
     * @brief Set gravity force
     *
     * @param x X component of gravity
     * @param y Y component of gravity
     */
    void setGravity(float x, float y) { gravity = {x, y}; }

    /**
     * @brief Set viscosity coefficient
     *
     * @param viscosity Viscosity coefficient value
     */
    void setViscosity(float viscosity) { viscosityCoefficient = viscosity; }

    /**
     * @brief Set gas constant for pressure calculation
     *
     * @param gasConstantValue Gas constant value
     */
    void setGasConstant(float gasConstantValue)
    {
      gasConstant = gasConstantValue;
    }

    /**
     * @brief Set rest density for pressure calculation
     *
     * @param density Rest density value
     */
    void setRestDensity(float density) { restDensity = density; }

    /**
     * @brief Set boundary damping factor
     *
     * @param damping Boundary damping coefficient
     */
    void setBoundaryDamping(float damping) { boundaryDamping = damping; }

    /**
     * @brief Compute density and pressure for each particle
     *
     * @param particleSystem Particle system with SoA layout
     * @param grid Spatial partitioning grid for neighbor search
     */
    void computeDensityPressure(ParticleSystem *particleSystem, Grid *grid);

    /**
     * @brief Compute forces (pressure, viscosity, gravity) for each particle
     *
     * @param particleSystem Particle system with SoA layout
     * @param grid Spatial partitioning grid for neighbor search
     */
    void computeForces(ParticleSystem *particleSystem, Grid *grid);

    /**
     * @brief Integrate particle positions and velocities
     *
     * @param particleSystem Particle system with SoA layout
     * @param dt Time step for integration
     */
    void integrate(ParticleSystem *particleSystem, float dt);

    /**
     * @brief Resolve collisions between particles and boundaries
     *
     * @param particleSystem Particle system with SoA layout
     * @param grid Spatial partitioning grid for neighbor search
     * @param width Width of simulation domain
     * @param height Height of simulation domain
     */
    void resolveCollisions(ParticleSystem *particleSystem, Grid *grid,
                           float width, float height);

    /**
     * @brief Compute forces exerted by boundaries on particles
     *
     * @param particleSystem Particle system with SoA layout
     * @param width Width of simulation domain
     * @param height Height of simulation domain
     */
    void computeBoundaryForces(ParticleSystem *particleSystem,
                               float width, float height);

    /**
     * @brief Compute density and pressure for a single particle by index
     *
     * @param particleSystem Particle system with SoA layout
     * @param particleIndex Index of the particle to process
     */
    void computeDensityPressureForParticle(ParticleSystem *particleSystem, size_t particleIndex);

    /**
     * @brief Compute forces (pressure, viscosity, gravity) for a single particle by index
     *
     * @param particleSystem Particle system with SoA layout
     * @param particleIndex Index of the particle to process
     */
    void computeForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex);

    /**
     * @brief Compute boundary forces for a single particle by index
     *
     * @param particleSystem Particle system with SoA layout
     * @param particleIndex Index of the particle to process
     * @param width Width of simulation domain
     * @param height Height of simulation domain
     */
    void computeBoundaryForcesForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height);

    /**
     * @brief Integrate position and velocity for a single particle by index
     *
     * @param particleSystem Particle system with SoA layout
     * @param particleIndex Index of the particle to process
     * @param dt Time step for integration
     */
    void integrateParticle(ParticleSystem *particleSystem, size_t particleIndex, float dt);

    /**
     * @brief Resolve collisions for a single particle by index
     *
     * @param particleSystem Particle system with SoA layout
     * @param particleIndex Index of the particle to process
     * @param width Width of simulation domain
     * @param height Height of simulation domain
     */
    void resolveCollisionsForParticle(ParticleSystem *particleSystem, size_t particleIndex, float width, float height);

    // Getter methods for SPH parameters
    float getSmoothingRadius() const { return h; }
    float getViscosity() const { return viscosityCoefficient; }
    float getGasConstant() const { return gasConstant; }
    float getRestDensity() const { return restDensity; }
    Vector2f getGravity() const { return gravity; }

    // Expose kernel functions for parallel tasks
    float kernelPoly6(float distSquared);
    Vector2f kernelGradSpiky(float distance, const Vector2f &direction);
    float kernelViscosityLaplacian(float distance);

  private:
    // SPH Parameters
    Vector2f gravity;                          // Gravity force
    float h;                                   // Smoothing radius
    float h2;                                  // h^2 (pre-computed)
    float viscosityCoefficient;                // Viscosity coefficient
    float gasConstant;                         // Gas constant for pressure
    float restDensity;                         // Rest density
    float boundaryDamping;                     // Boundary collision damping
    float gamma;                               // For equation of state
    int timeStepCounter;                       // For debugging
    std::unordered_set<size_t> debugParticles; // For tracking specific particles
  };

} // namespace sph
