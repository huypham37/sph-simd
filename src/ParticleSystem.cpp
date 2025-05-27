#include "ParticleSystem.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include "SPHConfig.hpp"

#ifndef HEADLESS_MODE
#include <SFML/Graphics.hpp>
#endif

namespace sph
{

	ParticleSystem::ParticleSystem(float width, float height, float smoothingRadius)
		: numParticles(0),
		  width(width),
		  height(height),
		  smoothingRadius(Config::SMOOTHING_RADIUS)
	{
		// Create Grid for spatial partitioning
		grid = new Grid(width, height, smoothingRadius);

		// Reserve initial capacity for better performance
		const size_t initialCapacity = 2500; // Slightly larger than typical particle count
		positionsX.reserve(initialCapacity);
		positionsY.reserve(initialCapacity);
		velocitiesX.reserve(initialCapacity);
		velocitiesY.reserve(initialCapacity);
		accelerationsX.reserve(initialCapacity);
		accelerationsY.reserve(initialCapacity);
		densities.reserve(initialCapacity);
		pressures.reserve(initialCapacity);
		masses.reserve(initialCapacity);
		cachedNeighbors.reserve(initialCapacity);

#ifndef HEADLESS_MODE
		shapes.reserve(initialCapacity);
		baseColors.reserve(initialCapacity);
#endif
	}

	ParticleSystem::~ParticleSystem()
	{
		reset();
		delete grid;
	}

	size_t ParticleSystem::addParticle(float x, float y)
	{
		// Add particle data to SoA
		positionsX.push_back(x);
		positionsY.push_back(y);
		velocitiesX.push_back(0.0f);
		velocitiesY.push_back(0.0f);
		accelerationsX.push_back(0.0f);
		accelerationsY.push_back(0.0f);
		densities.push_back(0.0f);
		pressures.push_back(0.0f);
		masses.push_back(Config::PARTICLE_MASS);
		cachedNeighbors.emplace_back(); // Empty neighbor list

#ifndef HEADLESS_MODE
		// Initialize visual components
		sf::CircleShape shape(Config::PARTICLE_RADIUS);
		sf::Color baseColor(40, 120, 255, 220);
		shape.setFillColor(baseColor);
		shape.setOrigin(sf::Vector2f(Config::PARTICLE_RADIUS, Config::PARTICLE_RADIUS));
		shape.setPosition(sf::Vector2f(x, y));
		shape.setOutlineThickness(0.5f);
		shape.setOutlineColor(sf::Color(10, 50, 150, 100));

		shapes.push_back(shape);
		baseColors.push_back(baseColor);
#endif

		return numParticles++;
	}

	void ParticleSystem::initializeDamBreak(int count)
	{
		// Clear any existing particles first
		reset();

		std::pair<int, int> ratio = Config::DAM_RATIO;
		float r = Config::PARTICLE_RADIUS;

		// Calculate dam dimensions in world units
		float damWidthUnits = ratio.second * 2 * r;
		float damHeightUnits = ratio.first * 2 * r;

		// Calculate number of particles per dimension to achieve desired count
		// while maintaining reasonable density
		float particleSpacing = Config::DP;
		int columns = static_cast<int>(damWidthUnits / particleSpacing);
		int rows = static_cast<int>(damHeightUnits / particleSpacing);

		// Adjust particle spacing to fit exact count if possible
		if (columns * rows > count)
		{
			// Reduce density if too many particles
			float ratio = std::sqrt(static_cast<float>(count) / (columns * rows));
			columns = static_cast<int>(columns * ratio);
			rows = static_cast<int>(rows * ratio);
			std::cout << "Adjusted particle grid to " << columns << "x" << rows << " = "
					  << (columns * rows) << " particles (requested " << count << ")" << std::endl;
		}

		float spacingX = damWidthUnits / columns;
		float spacingY = damHeightUnits / rows;

		// Small jitter for more natural appearance
		float jitter = particleSpacing * 0.1f;

		std::cout << "Creating dam break scenario with " << columns << "x" << rows << " particles" << std::endl;

		// Create the dam block of particles (positioned in the left side of the domain)
		for (int y = 0; y < rows; ++y)
		{
			for (int x = 0; x < columns; ++x)
			{
				float posX = x * spacingX;
				// Invert y-position to start from bottom
				float posY = height - (y + 1) * spacingY;

				// Ensure we don't exceed the desired particle count
				if (numParticles < static_cast<size_t>(count))
				{
					addParticle(posX, posY);
				}
				else
				{
					return;
				}
			}
		}

		std::cout << "Dam break scenario created with " << numParticles << " particles" << std::endl;
	}

	void ParticleSystem::reset()
	{
		// Clear all SoA vectors
		positionsX.clear();
		positionsY.clear();
		velocitiesX.clear();
		velocitiesY.clear();
		accelerationsX.clear();
		accelerationsY.clear();
		densities.clear();
		pressures.clear();
		masses.clear();
		cachedNeighbors.clear();

#ifndef HEADLESS_MODE
		shapes.clear();
		baseColors.clear();
#endif

		numParticles = 0;
	}

	void ParticleSystem::updateGrid()
	{
		grid->clear();
		// Update grid using direct SoA array access - much more cache-friendly
		const float *pos_x = positionsX.data();
		const float *pos_y = positionsY.data();

		// Process in batches for better cache locality
		constexpr size_t BATCH_SIZE = 64;
		for (size_t batch = 0; batch < numParticles; batch += BATCH_SIZE)
		{
			size_t batch_end = std::min(batch + BATCH_SIZE, numParticles);
			for (size_t i = batch; i < batch_end; ++i)
			{
				grid->insertParticleByIndex(i, pos_x[i], pos_y[i]);
			}
		}
	}

	const std::vector<size_t> &ParticleSystem::getCachedNeighbors(size_t index) const
	{
		return cachedNeighbors[index];
	}

	void ParticleSystem::setCachedNeighbors(size_t index, const std::vector<size_t> &neighbors)
	{
		cachedNeighbors[index] = neighbors;
	}

#ifndef HEADLESS_MODE
	void ParticleSystem::updateParticleVisuals(size_t index)
	{
		// Direct array access instead of expensive getters
		float vel_x = velocitiesX[index];
		float vel_y = velocitiesY[index];
		float speed = std::sqrt(vel_x * vel_x + vel_y * vel_y);

		// Color based on speed
		float normalizedSpeed = std::min(speed / 100.0f, 1.0f);
		std::uint8_t red = static_cast<std::uint8_t>(normalizedSpeed * 255);
		std::uint8_t blue = static_cast<std::uint8_t>((1.0f - normalizedSpeed) * 255);

		sf::Color currentColor(red, 120, blue, 220);
		shapes[index].setFillColor(currentColor);
		shapes[index].setPosition(sf::Vector2f(positionsX[index], positionsY[index]));
	}

	void ParticleSystem::updateAllParticleVisuals()
	{
		// Reduce visual update frequency - only update every 4th frame
		static int visualFrameCounter = 0;
		if (++visualFrameCounter % 4 != 0)
		{
			return; // Skip visual updates 3 out of 4 frames
		}

		// Use direct array access for better SIMD performance
		const float *vel_x = velocitiesX.data();
		const float *vel_y = velocitiesY.data();
		const float *pos_x = positionsX.data();
		const float *pos_y = positionsY.data();

		// Process in batches for better cache locality
		constexpr size_t BATCH_SIZE = 32; // Smaller batch for visual updates
		for (size_t batch = 0; batch < numParticles; batch += BATCH_SIZE)
		{
			size_t batch_end = std::min(batch + BATCH_SIZE, numParticles);

			// Vectorizable loop with direct array access
			for (size_t i = batch; i < batch_end; ++i)
			{
				float speed = std::sqrt(vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i]);
				float normalizedSpeed = std::min(speed / 100.0f, 1.0f);
				std::uint8_t red = static_cast<std::uint8_t>(normalizedSpeed * 255);
				std::uint8_t blue = static_cast<std::uint8_t>((1.0f - normalizedSpeed) * 255);

				sf::Color currentColor(red, 120, blue, 220);
				shapes[i].setFillColor(currentColor);
				shapes[i].setPosition(sf::Vector2f(pos_x[i], pos_y[i]));
			}
		}
	}
#endif

} // namespace sph