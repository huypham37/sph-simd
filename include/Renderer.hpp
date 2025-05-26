#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <memory>
#include "Particle.hpp"

namespace sph
{
	// Forward declaration
	class ParticleSystem;

	/**
	 * @brief Handles rendering of simulation elements
	 *
	 * Responsible for drawing particles.
	 */
	class Renderer
	{
	public:
		/**
		 * @brief Constructor
		 */
		Renderer();

		/**
		 * @brief Draw all particles
		 *
		 * @param particles Vector of particles to draw
		 * @param window SFML window to draw to
		 */
		void drawParticles(const std::vector<Particle *> &particles, sf::RenderWindow &window);

		/**
		 * @brief Draw particles from SoA particle system
		 *
		 * @param particleSystem Particle system using SoA layout
		 * @param window SFML window to draw to
		 */
		void drawParticles(const ParticleSystem *particleSystem, sf::RenderWindow &window);

	private:
		sf::Font font;	 // Font for text rendering
		bool fontLoaded; // Whether font was successfully loaded
	};

} // namespace sph