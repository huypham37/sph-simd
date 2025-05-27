#include "SPHSimulation.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <cmath> // Add cmath header for std::sqrt
#ifdef _OPENMP
#include <omp.h> // Include OpenMP header for thread functions
#endif

namespace sph
{

	SPHSimulation::SPHSimulation(float width, float height)
		: width(width),
		  height(height),
		  smoothingRadius(Config::SMOOTHING_RADIUS),
		  framesSinceLastMetricUpdate(0),
		  particlesPerSecond(0.0),
		  timeStepsPerSecond(0.0),
		  avgTimePerStep(0.0),
		  totalSimulationTime(0.0)
	{
		// Create component instances
		particles = std::make_unique<ParticleSystem>(width, height, smoothingRadius);
		physics = std::make_unique<SPHPhysics>();
		physics->setSmoothingRadius(smoothingRadius);
#ifndef HEADLESS_MODE
		renderer = std::make_unique<Renderer>();
#endif

		lastMetricUpdateTime = std::chrono::high_resolution_clock::now();

		std::cout << "SPH Simulation initialized with dimensions " << width << "x" << height << std::endl;
	}

	SPHSimulation::~SPHSimulation()
	{
		// Unique pointers will clean up automatically
	}

	void SPHSimulation::update(float dt)
	{
		// Start timing this step
		auto stepStartTime = std::chrono::high_resolution_clock::now();

		constexpr int sub_steps = 1; // Minimal substeps for stability
		float sub_dt = dt / sub_steps;

		// Update grid (needs to be sequential due to data structure updates)
		particles->updateGrid();

		// Substep loop for physics
		for (int step = 0; step < sub_steps; ++step)
		{
			size_t particleCount = particles->getParticleCount();

			// Batch processing with better cache locality
			constexpr size_t BATCH_SIZE = 64; // Optimized for cache line size

// Single OpenMP parallel region for the entire simulation step
#pragma omp parallel
			{
// Cache neighbors in parallel with batch processing
#pragma omp for schedule(dynamic, BATCH_SIZE / 4)
				for (size_t batch = 0; batch < particleCount; batch += BATCH_SIZE)
				{
					size_t batch_end = std::min(batch + BATCH_SIZE, particleCount);
					for (size_t i = batch; i < batch_end; ++i)
					{
						particles->setCachedNeighbors(i, particles->getGrid()->getNeighborsByIndex(i, smoothingRadius * 2, particles.get()));
					}
				}

// Density and pressure calculation phase with batch processing
#pragma omp barrier
#pragma omp for schedule(static, BATCH_SIZE)
				for (size_t i = 0; i < particleCount; ++i)
				{
					physics->computeDensityPressureForParticle(particles.get(), i);
				}

// Force calculation phase with batch processing
#pragma omp barrier
#pragma omp for schedule(static, BATCH_SIZE)
				for (size_t i = 0; i < particleCount; ++i)
				{
					physics->computeForcesForParticle(particles.get(), i);
				}

// Boundary forces with batch processing
#pragma omp barrier
#pragma omp for schedule(static, BATCH_SIZE)
				for (size_t i = 0; i < particleCount; ++i)
				{
					physics->computeBoundaryForcesForParticle(particles.get(), i, width, height);
				}

// Integration with batch processing
#pragma omp barrier
#pragma omp for schedule(static, BATCH_SIZE)
				for (size_t i = 0; i < particleCount; ++i)
				{
					physics->integrateParticle(particles.get(), i, sub_dt);
				}
			} // End of parallel region
		}

		// Collision resolution with batch processing for better cache locality
		size_t particleCount = particles->getParticleCount();
		constexpr size_t COLLISION_BATCH_SIZE = 128; // Larger batch for simple operations

#pragma omp parallel for schedule(static, COLLISION_BATCH_SIZE)
		for (size_t i = 0; i < particleCount; ++i)
		{
			physics->resolveCollisionsForParticle(particles.get(), i, width, height);
		}

		// Calculate metrics
		auto stepEndTime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> stepDuration = stepEndTime - stepStartTime;
		double stepTime = stepDuration.count();

		// Update metrics
		framesSinceLastMetricUpdate++;
		totalSimulationTime += stepTime;

		stepTimes.push_back(stepTime);
		if (stepTimes.size() > 100)
			stepTimes.erase(stepTimes.begin()); // Keep last 100 measurements

		// Calculate average time per step with better cache efficiency
		avgTimePerStep = 0.0;
		const double *step_times_data = stepTimes.data();
		size_t step_count = stepTimes.size();

// Vectorizable reduction with better memory access pattern
#pragma omp simd reduction(+ : avgTimePerStep)
		for (size_t i = 0; i < step_count; ++i)
		{
			avgTimePerStep += step_times_data[i];
		}
		avgTimePerStep /= step_count;

		// Update metrics every second
		auto now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = now - lastMetricUpdateTime;
		if (elapsed.count() >= 1.0)
		{ // Update metrics every second
			particlesPerSecond = particles->getParticleCount() * framesSinceLastMetricUpdate / elapsed.count();
			timeStepsPerSecond = framesSinceLastMetricUpdate / elapsed.count();

			framesSinceLastMetricUpdate = 0;
			lastMetricUpdateTime = now;
		}
	}

#ifndef HEADLESS_MODE
	void SPHSimulation::draw(sf::RenderWindow &window)
	{
		// Update visuals for all particles
		particles->updateAllParticleVisuals();

		// Draw particles using the new SoA interface
		renderer->drawParticles(particles.get(), window);
	}
#endif

	void SPHSimulation::reset()
	{
		particles->reset();
	}

	void SPHSimulation::initializeDefaultParticles(int count)
	{
		particles->initializeDamBreak(count);
	}

	size_t SPHSimulation::getParticleCount() const
	{
		return particles->getParticleCount();
	}

	void SPHSimulation::setGravity(float x, float y)
	{
		physics->setGravity(x, y);
	}

	void SPHSimulation::setViscosity(float v)
	{
		physics->setViscosity(v);
	}

	void SPHSimulation::setGasConstant(float k)
	{
		physics->setGasConstant(k);
	}

	void SPHSimulation::setRestDensity(float d)
	{
		physics->setRestDensity(d);
	}

	void SPHSimulation::applyMouseForce(const Vector2f &mousePos, float strength)
	{
		// Radius of influence for mouse interaction
		const float radiusOfInfluence = 50.0f;
		const float radiusSqr = radiusOfInfluence * radiusOfInfluence;

		size_t particleCount = particles->getParticleCount();

		// Direct array access for better cache performance
		std::vector<float> &vel_x = particles->getVelocitiesX();
		std::vector<float> &vel_y = particles->getVelocitiesY();
		const std::vector<float> &pos_x = particles->getPositionsX();
		const std::vector<float> &pos_y = particles->getPositionsY();

		// Process in batches for better cache locality
		constexpr size_t BATCH_SIZE = 64;
		for (size_t batch = 0; batch < particleCount; batch += BATCH_SIZE)
		{
			size_t batch_end = std::min(batch + BATCH_SIZE, particleCount);

// Vectorizable loop with direct array access
#pragma omp simd
			for (size_t i = batch; i < batch_end; ++i)
			{
				float dx = pos_x[i] - mousePos.x;
				float dy = pos_y[i] - mousePos.y;
				float distanceSqr = dx * dx + dy * dy;

				// Only apply force if particle is within radius of influence
				if (distanceSqr < radiusSqr)
				{
					// Calculate force based on distance (stronger when closer)
					float distance = std::sqrt(distanceSqr);

					// Normalize direction and scale by strength and distance factor
					if (distance > 0.1f) // Prevent division by zero or very small values
					{
						float inv_distance = 1.0f / distance;
						float dir_x = dx * inv_distance;
						float dir_y = dy * inv_distance;
						float forceFactor = strength * (1.0f - distance / radiusOfInfluence);

						// Apply force directly to velocity for immediate response
						vel_x[i] += dir_x * forceFactor;
						vel_y[i] += dir_y * forceFactor;
					}
				}
			}
		}
	}

	void SPHSimulation::printPerformanceMetrics() const
	{
		std::cout << "===== Performance Metrics =====" << std::endl;
		std::cout << "Particles: " << particles->getParticleCount() << std::endl;
#ifdef _OPENMP
		std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#else
		std::cout << "OpenMP: Not enabled" << std::endl;
#endif
		std::cout << "Particles/second: " << std::fixed << std::setprecision(1) << particlesPerSecond << std::endl;
		std::cout << "Timesteps/second: " << std::fixed << std::setprecision(2) << timeStepsPerSecond << std::endl;
		std::cout << "Timesteps/minute: " << std::fixed << std::setprecision(1) << timeStepsPerSecond * 60 << std::endl;
		std::cout << "Timesteps/hour: " << std::fixed << std::setprecision(1) << timeStepsPerSecond * 3600 << std::endl;
		std::cout << "Avg time per step: " << std::fixed << std::setprecision(5) << avgTimePerStep * 1000 << " ms" << std::endl;
		std::cout << "=============================" << std::endl;
	}

} // namespace sph