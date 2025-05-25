#pragma once

#include <memory>
#include <utility>

class SPHConfig
{

public:
  // Singleton Access
  static SPHConfig &getInstance()
  {
    static SPHConfig instance;
    return instance;
  };

  // Simulation space parameters
  static constexpr int WIDTH = 1200;
  static constexpr int HEIGHT = 600;
  static constexpr int PARTICLE_COUNT = 2000;
  static constexpr float TIME_STEP = 0.001f;                 // Important, please do not change this
  static constexpr std::pair<int, int> DAM_RATIO = {50, 40}; // this ratio following the particle count

  // SPH physics parameters
  static constexpr float DP = 10.0f;
  static constexpr float PARTICLE_RADIUS = DP * 0.5f;
  static constexpr float SMOOTHING_RADIUS = DP * 1.4f;
  static constexpr float GAS_CONSTANT = 5000000.0f;
  static constexpr float VISCOSITY = 0.1f;
  static constexpr float GRAVITY = 500.0f;
  static constexpr float REST_DENSITY = 0.9f;
  static constexpr float BOUNDARY_DAMPING = 0.3f;
  static constexpr float GAMMA = 7.0f;

  // Particle mass calculated as density * area (for 2D)
  static constexpr float PARTICLE_MASS = REST_DENSITY * (DP * DP);

  // Grid parameters
  static constexpr float CELL_SIZE = SMOOTHING_RADIUS; // Usually set to smoothingRadius

private:
  // private construtor for Singleton
  SPHConfig() {};

  // Delete copy construtor
  SPHConfig(const SPHConfig &) = delete;
  SPHConfig &operator=(const SPHConfig &) = delete;
};

// Using conviniently
using Config = SPHConfig;
