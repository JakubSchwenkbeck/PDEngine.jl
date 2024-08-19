# PDEngine.jl

**PDEngine.jl** is a Julia library designed for solving partial differential equations (PDEs) using a variety of numerical methods. This library aims to provide efficient, flexible, and easy-to-use solvers for some of the most common PDEs encountered in scientific and engineering problems.

## Features

- **Heat Equation**: Models heat distribution over time.
- **Wave Equation**: Simulates the propagation of waves through a medium.
- **Navier-Stokes Equations**: Governs the motion of viscous fluid substances.
- **Poisson's Equation**: Used in electrostatics, mechanical engineering, and theoretical physics.

## Numerical Methods

- **Finite Difference Method (FDM)**: A straightforward approach to discretize PDEs on regular grids.
- **Finite Element Method (FEM)**: Ideal for complex geometries and adaptive mesh refinement.
- **Spectral Methods**: Provides high accuracy for smooth solutions and periodic boundary conditions.

## Getting Started

To use PDEngine.jl, you need to have Julia installed. You can install the package using Julia's package manager.

1. **Install PDEngine.jl**:
    ```julia
    using Pkg
    Pkg.add("PDEngine")
    ```

2. **Using the Package**:
    ```julia
    using PDEngine

    # Example: Solving the heat equation in 1D
    Δx = 0.01
    Δt = 0.0001
    T = 0.1
    α = 0.01
    N = 100
    temperature_distribution = heat_eq_1d(N, α, T, Δx, Δt)

    # Example: Solving Poisson's equation in 2D
    f = zeros(N+1, N+1)  # Example source term
    Δx = 0.01
    poisson_solution = poisson_2d(N, f, Δx)
    ```
