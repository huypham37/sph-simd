#include "Particle.hpp"
#include <cmath>
#include "SPHConfig.hpp"

namespace sph {

Particle::Particle(float x, float y)
	: position({x, y}),
	  velocity({0.0f, 0.0f}),
	  acceleration({0.0f, 0.0f}),
	  mass(Config::PARTICLE_MASS),
	  density(0.0f),
	  pressure(0.0f)
#ifndef HEADLESS_MODE
	  ,shape(Config::PARTICLE_RADIUS)
#endif
{
#ifndef HEADLESS_MODE
	// Base color for still particles - only in graphical mode
	baseColor = sf::Color(40, 120, 255, 220);
	shape.setFillColor(baseColor);
	shape.setOrigin({RADIUS, RADIUS});
	shape.setPosition(position);

	// Add slight outline for better visibility
	shape.setOutlineThickness(0.5f);
	shape.setOutlineColor(sf::Color(10, 50, 150, 100));
#endif
}

void Particle::update(float dt)
{
	// Update position
	position += velocity * dt;

#ifndef HEADLESS_MODE
	// Update color based on velocity - but only every few frames
	static int frameCount = 0;
	frameCount = (frameCount + 1) % 20; // Update color every 5 frames
	if (frameCount == 0)
	{
		updateColor();
	}

	// Update shape position for rendering
	shape.setPosition(position);
#endif
}

#ifndef HEADLESS_MODE
void Particle::updateVisuals()
{
	updateColor();
}

void Particle::draw(sf::RenderWindow &window)
{
	window.draw(shape);
}

void Particle::updateColor()
{
	// Use squared speed to avoid sqrt when possible
	float speedSq = velocity.x * velocity.x + velocity.y * velocity.y;

	// Early exit for very slow particles
	if (speedSq < 1.0f)
	{
		// If barely moving, just use base color (avoid expensive calculations)
		shape.setFillColor(baseColor);
		shape.setOutlineColor(sf::Color(10, 50, 150, 100));
		shape.setOutlineThickness(0.5f);

		// Reset radius if not at default
		float currentRadius = shape.getRadius();
		if (std::abs(currentRadius - RADIUS) > 0.01f)
		{
			shape.setRadius(RADIUS);
			shape.setOrigin({RADIUS, RADIUS});
		}
		return;
	}

	// For moving particles, compute the full color gradient
	float speed = std::sqrt(speedSq);
	constexpr float MAX_SPEED = 100.0f;
	float speedFactor = speed / MAX_SPEED;
	if (speedFactor > 1.0f)
		speedFactor = 1.0f;

	// Map speed to a color gradient
	// - Low speed: blue (cold)
	// - Medium speed: light blue/cyan
	// - High speed: white with blue tint

	int r = baseColor.r + static_cast<int>((255 - baseColor.r) * speedFactor * 0.8f);
	int g = baseColor.g + static_cast<int>((255 - baseColor.g) * speedFactor);
	int b = baseColor.b + static_cast<int>((255 - baseColor.b) * speedFactor * 0.4f);

	// Slight glow effect for fast-moving particles
	int alpha = baseColor.a + static_cast<int>((255 - baseColor.a) * speedFactor * 0.3f);

	shape.setFillColor(sf::Color(r, g, b, alpha));

	// Adjust radius slightly based on speed for more dynamic appearance
	float radiusMultiplier = 1.0f + speedFactor * 0.15f;
	shape.setRadius(RADIUS * radiusMultiplier);
	shape.setOrigin({RADIUS * radiusMultiplier, RADIUS * radiusMultiplier});

	// Update outline based on speed as well
	if (speedFactor > 0.7f)
	{
		shape.setOutlineColor(sf::Color(200, 230, 255, 120));
		shape.setOutlineThickness(0.8f);
	}
	else
	{
		shape.setOutlineColor(sf::Color(10, 50, 150, 100));
		shape.setOutlineThickness(0.5f);
	}
}
#endif // HEADLESS_MODE

} // namespace sph