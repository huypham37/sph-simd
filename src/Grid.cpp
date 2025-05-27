#include "Grid.hpp"
#include "ParticleSystem.hpp"
#include <cmath>

namespace sph
{

  Grid::Grid(float width, float height, float cellSize)
      : cellSize(cellSize), invCellSize(1.0f / cellSize),
        gridWidth(static_cast<int>(std::ceil(width / cellSize))),
        gridHeight(static_cast<int>(std::ceil(height / cellSize))) {}

  void Grid::clear() { cells.clear(); }

  void Grid::insertParticleByIndex(size_t particleIndex, float x, float y)
  {
    int index = getCellIndex(x, y);
    cells[index].push_back(particleIndex);
  }

  void Grid::updateGrid(const ParticleSystem *particleSystem)
  {
    this->clear();

    // Use direct array access instead of getPosition for better performance
    const float *pos_x = particleSystem->getPositionsXData();
    const float *pos_y = particleSystem->getPositionsYData();

    for (size_t i = 0; i < particleSystem->getParticleCount(); ++i)
    {
      insertParticleByIndex(i, pos_x[i], pos_y[i]);
    }
  }

  std::vector<size_t> Grid::getNeighbors(float x, float y, float radius, const ParticleSystem *particleSystem)
  {
    std::vector<size_t> neighbors;

    // Calculate cell coordinates and radius in cells
    auto [cellX, cellY] = positionToCell(x, y);
    int cellRadius = static_cast<int>(std::ceil(radius * invCellSize));
    float radiusSquared = radius * radius; // For distance check optimization

    // Iterate through potential neighbor cells
    for (int i = -cellRadius; i <= cellRadius; ++i)
    {
      for (int j = -cellRadius; j <= cellRadius; ++j)
      {
        int neighborCellX = cellX + i;
        int neighborCellY = cellY + j;

        // Skip if outside grid bounds
        if (neighborCellX < 0 || neighborCellX >= gridWidth ||
            neighborCellY < 0 || neighborCellY >= gridHeight)
        {
          continue;
        }

        int cellIdx = getCellIndex(neighborCellX, neighborCellY);
        auto it = cells.find(cellIdx);
        if (it != cells.end())
        {
          // Filter particles by actual distance - use direct array access
          const float *pos_x = particleSystem->getPositionsXData();
          const float *pos_y = particleSystem->getPositionsYData();

          for (size_t particleIndex : it->second)
          {
            float dx = x - pos_x[particleIndex];
            float dy = y - pos_y[particleIndex];
            float distSquared = dx * dx + dy * dy;

            if (distSquared <= radiusSquared)
            {
              neighbors.push_back(particleIndex);
            }
          }
        }
      }
    }

    return neighbors;
  }

  std::vector<size_t> Grid::getNeighborsByIndex(size_t particleIndex, float radius, const ParticleSystem *particleSystem)
  {
    // Use direct array access instead of getPosition
    const float *pos_x = particleSystem->getPositionsXData();
    const float *pos_y = particleSystem->getPositionsYData();
    return getNeighbors(pos_x[particleIndex], pos_y[particleIndex], radius, particleSystem);
  }

  int Grid::getCellIndex(float x, float y) const
  {
    auto [cellX, cellY] = positionToCell(x, y);
    return getCellIndex(cellX, cellY);
  }

  int Grid::getCellIndex(int cellX, int cellY) const
  {
    return cellY * gridWidth + cellX;
  }

  std::pair<int, int> Grid::positionToCell(float x, float y) const
  {
    int cellX = static_cast<int>(x * invCellSize);
    int cellY = static_cast<int>(y * invCellSize);
    return {cellX, cellY};
  }

} // namespace sph
