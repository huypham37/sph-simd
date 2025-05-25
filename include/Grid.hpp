#pragma once

#include <vector>
#include <unordered_map>

#ifdef HEADLESS_MODE
#include <SFML/System.hpp> // Use the main System header, not specific files
#else
#include <SFML/Graphics.hpp>
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