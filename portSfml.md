# SPH Fluid Simulation - SFML 3.0 Migration Guide

This document outlines the changes made to port the SPH Fluid Simulation from older SFML versions to SFML 3.0.

## Key Changes for SFML 3.0 Compatibility

### CMake Configuration
- Updated target linking from `sfml-graphics` to `SFML::Graphics`
- Updated SFML version requirement from 2.5 to 3.0

### Vector Initialization
- Changed from constructor style to braced initialization:
  - Old: `sf::Vector2f(x, y)`
  - New: `sf::Vector2f{x, y}` or `{x, y}`

### Text Initialization
- Updated the text constructor:
  - Old: `sf::Text text; text.setFont(font);` 
  - New: `sf::Text text(font, "");`

### VideoMode Initialization
- Updated VideoMode initialization:
  - Old: `sf::VideoMode(width, height)`
  - New: `sf::VideoMode({width, height})`

### Event Handling
- Changed event handling to use std::optional and type-safe event access:
  - Old: 
    ```cpp
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed) {
            // ...
        }
    }
    ```
  - New:
    ```cpp
    std::optional<sf::Event> eventOpt;
    while ((eventOpt = window.pollEvent())) {
        const auto& event = *eventOpt;
        if (event.is<sf::Event::Closed>()) {
            // ...
        }
    }
    ```

- Updated event member access:
  - Old: `event.mouseButton.button == sf::Mouse::Left`
  - New: `const auto* mouseEvent = event.getIf<sf::Event::MouseButtonPressed>(); mouseEvent->button == sf::Mouse::Button::Left`

### Enum Classes
- Updated enums to use scoped enum syntax:
  - Old: `sf::Keyboard::R`
  - New: `sf::Keyboard::Key::R`
  - Old: `sf::Mouse::Left`
  - New: `sf::Mouse::Button::Left`

### Font Loading
- Changed font loading function:
  - Old: `font.loadFromFile()`
  - New: `font.openFromFile()`

## Additional Notes
- SFML 3.0 has improved type safety and better compile-time errors
- The API is more consistent with modern C++ practices
- Many performance improvements in rendering and event handling
- Enhanced support for high DPI displays

## Building the Project
To build the SPH Fluid Simulation with SFML 3.0:

1. Make sure SFML 3.0 is installed on your system
2. Run the build script: `./build.sh`
3. Execute the simulation: `./build/sph_simulation`