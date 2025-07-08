# Quantum Proton Simulation

![Proton Simulation Demo](demo.gif) *(Include an actual gif/screenshot after uploading)*

A Python/Pygame visualization of quark confinement in protons, simulating the strong nuclear force's unique behavior:
- Short-range repulsion
- Intermediate attraction
- Long-range fading (confinement)

## Features
- **Dynamic spring system** connecting quarks with realistic physics
- **Modified Yukawa potential** modeling the strong force
- **Interactive controls** to inject energy into specific quarks
- **Visual effects** including particle trails, glow effects, and tension-sensitive springs
- **Educational UI** explaining the physics model

## Physics Model
Uses a hybrid approach combining:
```python
F(r) = A/r² - B·exp(-r/D)/r - C·exp(-r/0.5)
