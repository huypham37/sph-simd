cmake_minimum_required(VERSION 3.14)
project(sph_simulation VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# OpenMP configuration for macOS/Apple Clang
if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "Clang")
    # Check if OpenMP is installed via Homebrew
    if(EXISTS "/opt/homebrew/opt/libomp")
        # OpenMP installed via Homebrew on Apple Silicon Mac
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    elseif(EXISTS "/usr/local/opt/libomp")
        # OpenMP installed via Homebrew on Intel Mac
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    else()
        message(WARNING "OpenMP not found via Homebrew. Please install with: brew install libomp")
    endif()
endif()

# Option for headless mode
option(HEADLESS_MODE "Build for headless mode without SFML graphics" OFF)

# Common source files
set(SOURCES
    src/main.cpp
    src/SPHSimulation.cpp
    src/SPHPhysics.cpp
    src/Particle.cpp
    src/Grid.cpp
    src/ParticleSystem.cpp
)

# Add Renderer.cpp only when not in headless mode
if(NOT HEADLESS_MODE)
    list(APPEND SOURCES src/Renderer.cpp)
endif()

# Create executable
add_executable(sph_simulation ${SOURCES})
target_include_directories(sph_simulation PRIVATE include)

# Add SIMD optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Clang-specific optimization flags
    target_compile_options(sph_simulation PRIVATE 
        -march=native          # Use all available CPU instructions
        -mtune=native          # Optimize for current CPU
        -ffast-math            # Enable fast math optimizations
        -funroll-loops         # Unroll loops for vectorization
        -fvectorize            # Enable auto-vectorization (Clang)
        -ftree-vectorize       # Enable auto-vectorization
    )
    
    # Add SIMD-specific flags for Apple Silicon (M1/M2)
    if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        target_compile_options(sph_simulation PRIVATE -mcpu=apple-m1)
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # GCC-specific optimization flags
    target_compile_options(sph_simulation PRIVATE 
        -march=native          # Use all available CPU instructions
        -mtune=native          # Optimize for current CPU
        -ffast-math            # Enable fast math optimizations
        -funroll-loops         # Unroll loops for vectorization
        -ftree-vectorize       # Enable auto-vectorization (GCC)
        -ftree-loop-vectorize  # Enable loop vectorization (GCC)
    )
endif()

# macOS specific SFML paths
if(APPLE)
    # Apple Silicon Mac
    if(EXISTS "/opt/homebrew/include")
        target_include_directories(sph_simulation PRIVATE /opt/homebrew/include)
        link_directories(/opt/homebrew/lib)
    endif()
    
    # Intel Mac
    if(EXISTS "/usr/local/include")
        target_include_directories(sph_simulation PRIVATE /usr/local/include)
        link_directories(/usr/local/lib)
    endif()
endif()

# Find SFML packages
if(HEADLESS_MODE)
    # Define HEADLESS_MODE for conditional compilation
    target_compile_definitions(sph_simulation PRIVATE HEADLESS_MODE)
    
    # Only find and link against SFML System for headless mode
    find_package(SFML 3.0 COMPONENTS System REQUIRED)

    # Add explicit library and linking info
    if(APPLE)
        if(EXISTS "/opt/homebrew/lib/libsfml-system.dylib")
            target_link_libraries(sph_simulation PRIVATE "/opt/homebrew/lib/libsfml-system.dylib")
        elseif(EXISTS "/usr/local/lib/libsfml-system.dylib")
            target_link_libraries(sph_simulation PRIVATE "/usr/local/lib/libsfml-system.dylib")
        else()
            target_link_libraries(sph_simulation PRIVATE sfml-system)
        endif()
    else()
        target_link_libraries(sph_simulation PRIVATE sfml-system)
    endif()
else()
    # Find and link against full SFML for graphical mode
    find_package(SFML 3.0 COMPONENTS System Window Graphics REQUIRED)

    # Add explicit library paths for Apple systems
    if(APPLE)
        if(EXISTS "/opt/homebrew/lib/libsfml-system.dylib" AND 
           EXISTS "/opt/homebrew/lib/libsfml-window.dylib" AND
           EXISTS "/opt/homebrew/lib/libsfml-graphics.dylib")
            target_link_libraries(sph_simulation PRIVATE 
                "/opt/homebrew/lib/libsfml-system.dylib"
                "/opt/homebrew/lib/libsfml-window.dylib"
                "/opt/homebrew/lib/libsfml-graphics.dylib"
            )
        elseif(EXISTS "/usr/local/lib/libsfml-system.dylib" AND 
               EXISTS "/usr/local/lib/libsfml-window.dylib" AND
               EXISTS "/usr/local/lib/libsfml-graphics.dylib")
            target_link_libraries(sph_simulation PRIVATE 
                "/usr/local/lib/libsfml-system.dylib"
                "/usr/local/lib/libsfml-window.dylib"
                "/usr/local/lib/libsfml-graphics.dylib"
            )
        else()
            target_link_libraries(sph_simulation PRIVATE sfml-system sfml-window sfml-graphics)
        endif()
    else()
        target_link_libraries(sph_simulation PRIVATE sfml-system sfml-window sfml-graphics)
    endif()
endif()

# OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sph_simulation PRIVATE OpenMP::OpenMP_CXX)
endif()

# Print some useful information during configure
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Headless mode: ${HEADLESS_MODE}")
message(STATUS "SFML include dirs: ${SFML_INCLUDE_DIR}")

# Install target
install(TARGETS sph_simulation DESTINATION bin)