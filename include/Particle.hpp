#pragma once

#ifdef HEADLESS_MODE
    #include <SFML/System.hpp>  // Use the main System header, not specific files
#else
    #include <SFML/Graphics.hpp>
#endif

#include <vector>

namespace sph {

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

	sf::Vector2f getPosition() const { return position; }
	sf::Vector2f getVelocity() const { return velocity; }
	sf::Vector2f getImmediateVelocity() const { return immediateVel; }
	sf::Vector2f getAcceleration() const { return acceleration; }
	float getMass() const { return mass; }
	float getDensity() const { return density; }
	float getPressure() const { return pressure; }

	void setPosition(const sf::Vector2f &pos) { position = pos; }
	void setVelocity(const sf::Vector2f &vel) { velocity = vel; }
	void setAcceleration(const sf::Vector2f &acc) { acceleration = acc; }
	void setDensity(float d) { density = d; }
	void setPressure(float p) { pressure = p; }
	void setImmediateVelocity(const sf::Vector2f &vel) {
		immediateVel = vel;}
	

	private:
		sf::Vector2f position;
		sf::Vector2f velocity;
		sf::Vector2f immediateVel;
		sf::Vector2f acceleration;
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