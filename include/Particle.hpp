#pragma once

#include "Grid.hpp" // For Vector2f abstraction

#ifndef HEADLESS_MODE
#include <SFML/Graphics.hpp>
#endif

#include <vector>

namespace sph
{

	class Particle
	{
	public:
		Particle(float x, float y);
		~Particle() = default;

		std::vector<Particle *> cachedNeighbors;
		void update(float dt);

#ifndef HEADLESS_MODE
		// Update the particle color based on velocity - only in graphical mode
		void updateColor();
		void updateVisuals();
		void draw(sf::RenderWindow &window);
#endif

		Vector2f getPosition() const { return position; }
		Vector2f getVelocity() const { return velocity; }
		Vector2f getImmediateVelocity() const { return immediateVel; }
		Vector2f getAcceleration() const { return acceleration; }
		float getMass() const { return mass; }
		float getDensity() const { return density; }
		float getPressure() const { return pressure; }

		void setPosition(const Vector2f &pos) { position = pos; }
		void setVelocity(const Vector2f &vel) { velocity = vel; }
		void setAcceleration(const Vector2f &acc) { acceleration = acc; }
		void setDensity(float d) { density = d; }
		void setPressure(float p) { pressure = p; }
		void setImmediateVelocity(const Vector2f &vel)
		{
			immediateVel = vel;
		}

	private:
		Vector2f position;
		Vector2f velocity;
		Vector2f immediateVel;
		Vector2f acceleration;
		float mass;
		float density;
		float pressure;

#ifndef HEADLESS_MODE
		sf::CircleShape shape;
		sf::Color baseColor; // Base color for the particle
#endif

		static constexpr float RADIUS = 5.0f;
	};

} // namespace sph