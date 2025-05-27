#pragma once

#ifndef HEADLESS_MODE
#include <SFML/Graphics.hpp>
#include "Renderer.hpp"
#endif

#include <memory>
#include <chrono>
#include <vector>
#include "ParticleSystem.hpp"
#include "SPHPhysics.hpp"

namespace sph
{

	/**
	 * @brief Main coordinator class for the SPH simulation
	 *
	 * Manages the overall simulation, delegating specific responsibilities
	 * to specialized components like physics, particle management,
	 * and rendering.
	 */
	class SPHSimulation
	{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param width Width of simulation domain
		 * @param height Height of simulation domain
		 */
		SPHSimulation(float width, float height);

		/**
		 * @brief Destructor
		 */
		~SPHSimulation();

		/**
		 * @brief Update simulation state for one time step
		 *
		 * @param dt Time step
		 */
		void update(float dt);

#ifndef HEADLESS_MODE
		/**
		 * @brief Render simulation to SFML window
		 *
		 * @param window SFML window to draw to
		 */
		void draw(sf::RenderWindow &window);
#endif

		/**
		 * @brief Initialize simulation with default particle configuration
		 *
		 * @param count Number of particles to create
		 */
		void initializeDefaultParticles(int count);

		/**
		 * @brief Add single particle at specific position
		 *
		 * @param x X-coordinate
		 * @param y Y-coordinate
		 */
		void addParticle(float x, float y);

		/**
		 * @brief Add multiple particles
		 *
		 * @param count Number of particles to add
		 */
		void addParticles(int count);

		/**
		 * @brief Remove specific number of particles
		 *
		 * @param count Number of particles to remove
		 */
		void removeParticles(int count);

		/**
		 * @brief Reset simulation (remove all particles)
		 */
		void reset();

		/**
		 * @brief Apply force to particles near mouse position
		 *
		 * @param mousePos Mouse position in window coordinates
		 * @param strength Force strength (positive = push, negative = pull)
		 */
		void applyMouseForce(const Vector2f &mousePos, float strength);

		// Simulation parameter setters
		void setGravity(float x, float y);
		void setViscosity(float v);
		void setGasConstant(float k);
		void setRestDensity(float d);

		// Statistics
		size_t getParticleCount() const;

		// New methods for metrics
		void printPerformanceMetrics() const;

		// Metrics getters
		double getParticlesPerSecond() const { return particlesPerSecond; }
		double getTimeStepsPerSecond() const { return timeStepsPerSecond; }
		double getAvgTimePerStep() const { return avgTimePerStep; }

	private:
		// Main simulation components
		std::unique_ptr<ParticleSystem> particles;
		std::unique_ptr<SPHPhysics> physics;
		std::unique_ptr<Grid> grid;
#ifndef HEADLESS_MODE
		std::unique_ptr<Renderer> renderer;
#endif

		// Simulation dimensions
		float width;
		float height;

		// SPH parameters
		float smoothingRadius;

		// Performance metrics
		std::chrono::time_point<std::chrono::high_resolution_clock> lastMetricUpdateTime;
		int framesSinceLastMetricUpdate;
		double particlesPerSecond;
		double timeStepsPerSecond;
		double avgTimePerStep;
		double totalSimulationTime;
		std::vector<double> stepTimes; // For tracking step time vs particle count
	};

} // namespace sph