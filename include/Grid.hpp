#pragma once

#include <vector>
#include <unordered_map>
#include <cstddef> // For size_t

#ifdef HEADLESS_MODE
// Simple Vector2f replacement for headless mode
struct Vector2f {
    float x, y;
    Vector2f(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}
    
    Vector2f operator-(const Vector2f& other) const {
        return Vector2f(x - other.x, y - other.y);
    }
    
    Vector2f operator+(const Vector2f& other) const {
        return Vector2f(x + other.x, y + other.y);
    }
    
    Vector2f operator*(float scalar) const {
        return Vector2f(x * scalar, y * scalar);
    }
    
    Vector2f operator/(float scalar) const {
        return Vector2f(x / scalar, y / scalar);
    }
    
    Vector2f& operator+=(const Vector2f& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    
    Vector2f& operator/=(float scalar) {
        x /= scalar;
        y /= scalar;
        return *this;
    }
};

// Global operator for scalar * Vector2f  
inline Vector2f operator*(float scalar, const Vector2f& vec) {
    return Vector2f(vec.x * scalar, vec.y * scalar);
}
#else
#include <SFML/Graphics.hpp>
using Vector2f = sf::Vector2f;
#endif

namespace sph
{
	// Forward declaration
	class ParticleSystem;

	class Grid
	{
	public:
		Grid(float width, float height, float cellSize);
		~Grid() = default;

		void clear();
		void insertParticleByIndex(size_t particleIndex, float x, float y);
		std::vector<size_t> getNeighbors(float x, float y, float radius, const ParticleSystem *particleSystem);
		std::vector<size_t> getNeighborsByIndex(size_t particleIndex, float radius, const ParticleSystem *particleSystem);
		void updateGrid(const ParticleSystem *particleSystem);
		// Add a getter for cell size to make its use explicit
		float getCellSize() const { return cellSize; }

	private:
		float cellSize;
		float invCellSize; // Inverse of cell size for optimization
		int gridWidth;
		int gridHeight;

		// Grid cell to particle indices mapping
		std::unordered_map<int, std::vector<size_t>> cells;

		// Helper methods
		int getCellIndex(float x, float y) const;
		int getCellIndex(int cellX, int cellY) const;
		std::pair<int, int> positionToCell(float x, float y) const;
	};

} // namespace sph