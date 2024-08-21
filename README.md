# PDEngine.jl (v0.0.9)

**PDEngine.jl** is a Julia library designed for solving partial differential equations (PDEs) using a variety of numerical methods. This library aims to provide efficient, flexible, and easy-to-use solvers for some of the most common PDEs encountered in scientific and engineering problems.
[Documentation](https://jakubschwenkbeck.github.io/PDEngine/)

## v0.0.9 now available : "pkg> add PDEngine" !!

## Features

- **Heat Equation**: Models heat distribution over time.
- **Wave Equation**: Simulates the propagation of waves through a medium.
- **Navier-Stokes Equations**: Governs the motion of viscous fluid substances.
- **Poisson's Equation**: Used in electrostatics, mechanical engineering, and theoretical physics.

## Numerical Methods

- **Finite Difference Method (FDM)**: A straightforward approach to discretize PDEs on regular grids.
- **Finite Element Method (FEM)**: Ideal for complex geometries and adaptive mesh refinement.
- **Spectral Methods**: Provides high accuracy for smooth solutions and periodic boundary conditions.

## Usage
**Using the Package**:
```julia
    using PDEngine

    # Example: Solving the heat equation in 1D
    Δx = 0.01
    Δt = 0.0001
    T = 0.1
    α = 0.01
    N = 100
    temperature_distribution = heat_eq_1d_fdm(N, α, T, Δx, Δt) # solving the heat equation with the finite differences method

    # Example: Solving Poisson's equation in 2D
    f = zeros(N+1, N+1)  # Example source term
    Δx = 0.01
    poisson_solution = poisson_2d_fem(N, f, Δx)  # solving the poisson equations with the finite elements method
```

## Currently working on:
### v0.0.95
- Better Function Naming
- Impementing Spectral methods
- Setting a default method with very simple Name like heat(params)
### v0.1.0
- Parallel computing to get the full speed of julia
- more methods and or PEDs
- ensure all Methods work without any flaws
- ...
    
