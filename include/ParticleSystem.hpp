#pragma once

#include <vector>
#include "Grid.hpp"

#ifndef HEADLESS_MODE
#include <SFML/Graphics.hpp>
#endif

namespace sph
{

	/**
	 * @brief Manages system of particles using Structure of Arrays (SoA)
	 *
	 * Responsible for managing all particles in the simulation using SoA layout
	 * for improved SIMD performance, their creation, deletion, and spatial partitioning with Grid.
	 */
	class ParticleSystem
	{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param width Width of simulation domain
		 * @param height Height of simulation domain
		 * @param smoothingRadius SPH smoothing radius
		 */
		ParticleSystem(float width, float height, float smoothingRadius);

		/**
		 * @brief Destructor
		 */
		~ParticleSystem();

		/**
		 * @brief Add particle at specific position
		 *
		 * @param x X-coordinate
		 * @param y Y-coordinate
		 * @return Index of the created particle
		 */
		size_t addParticle(float x, float y);

		/**
		 * @brief Initialize particles in a dam break configuration
		 *
		 * Creates a dense block of particles in one side of the domain
		 * that will collapse and flow when simulation starts
		 *
		 * @param count Number of particles to create
		 */
		void initializeDamBreak(int count);

		/**
		 * @brief Delete all particles
		 */
		void reset();

		/**
		 * @brief Update grid for spatial partitioning
		 */
		void updateGrid();

		/**
		 * @brief Get total number of particles
		 *
		 * @return Particle count
		 */
		size_t getParticleCount() const { return numParticles; }

		/**
		 * @brief Get simulation domain dimensions
		 *
		 * @return Width of simulation domain
		 */
		float getWidth() const { return width; }

		/**
		 * @brief Get simulation domain dimensions
		 *
		 * @return Height of simulation domain
		 */
		float getHeight() const { return height; }

		/**
		 * @brief Get pointer to spatial partitioning grid
		 *
		 * @return Pointer to Grid
		 */
		Grid *getGrid() const { return grid; }

		// SoA accessors - for accessing individual particle data arrays
		const std::vector<float> &getPositionsX() const { return positionsX; }
		const std::vector<float> &getPositionsY() const { return positionsY; }
		const std::vector<float> &getVelocitiesX() const { return velocitiesX; }
		const std::vector<float> &getVelocitiesY() const { return velocitiesY; }
		const std::vector<float> &getAccelerationsX() const { return accelerationsX; }
		const std::vector<float> &getAccelerationsY() const { return accelerationsY; }
		const std::vector<float> &getDensities() const { return densities; }
		const std::vector<float> &getPressures() const { return pressures; }
		const std::vector<float> &getMasses() const { return masses; }

		// Non-const accessors for physics calculations
		std::vector<float> &getPositionsX() { return positionsX; }
		std::vector<float> &getPositionsY() { return positionsY; }
		std::vector<float> &getVelocitiesX() { return velocitiesX; }
		std::vector<float> &getVelocitiesY() { return velocitiesY; }
		std::vector<float> &getAccelerationsX() { return accelerationsX; }
		std::vector<float> &getAccelerationsY() { return accelerationsY; }
		std::vector<float> &getDensities() { return densities; }
		std::vector<float> &getPressures() { return pressures; }
		std::vector<float> &getMasses() { return masses; }

		// Helper methods for getting particle data by index
		Vector2f getPosition(size_t index) const;
		Vector2f getVelocity(size_t index) const;
		Vector2f getAcceleration(size_t index) const;
		float getDensity(size_t index) const;
		float getPressure(size_t index) const;
		float getMass(size_t index) const;

		void setPosition(size_t index, const Vector2f &pos);
		void setVelocity(size_t index, const Vector2f &vel);
		void setAcceleration(size_t index, const Vector2f &acc);
		void setDensity(size_t index, float density);
		void setPressure(size_t index, float pressure);

		// Cached neighbors accessors
		const std::vector<size_t> &getCachedNeighbors(size_t index) const;
		void setCachedNeighbors(size_t index, const std::vector<size_t> &neighbors);

#ifndef HEADLESS_MODE
		const std::vector<sf::CircleShape> &getShapes() const { return shapes; }
		std::vector<sf::CircleShape> &getShapes() { return shapes; }
		const std::vector<sf::Color> &getBaseColors() const { return baseColors; }
		std::vector<sf::Color> &getBaseColors() { return baseColors; }

		void updateParticleVisuals(size_t index);
		void updateAllParticleVisuals();
#endif

	private:
		// Structure of Arrays (SoA) for particle data
		std::vector<float> positionsX;
		std::vector<float> positionsY;
		std::vector<float> velocitiesX;
		std::vector<float> velocitiesY;
		std::vector<float> accelerationsX;
		std::vector<float> accelerationsY;
		std::vector<float> densities;
		std::vector<float> pressures;
		std::vector<float> masses;

		// Neighbor cache using indices instead of pointers
		std::vector<std::vector<size_t>> cachedNeighbors;

#ifndef HEADLESS_MODE
		// Visual data (only for graphical mode)
		std::vector<sf::CircleShape> shapes;
		std::vector<sf::Color> baseColors;
#endif

		size_t numParticles;   // Current number of particles
		Grid *grid;			   // Spatial partitioning grid
		float width;		   // Width of simulation domain
		float height;		   // Height of simulation domain
		float smoothingRadius; // SPH smoothing radius
	};

} // namespace sph