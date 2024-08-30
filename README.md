# PDEngine.jl [![version](https://juliahub.com/docs/General/PDEngine/stable/version.svg)](https://juliahub.com/ui/Packages/General/PDEngine)
**PDEngine.jl** is a Julia library designed for solving partial differential equations (PDEs) using a variety of numerical methods. This library aims to provide efficient, flexible, and easy-to-use solvers for some of the most common PDEs encountered in scientific and engineering problems.
[Documentation](https://jakubschwenkbeck.github.io/PDEngine.jl/)

## Features

- **Heat Equation**: Models heat distribution over time.
- **Wave Equation**: Simulates the propagation of waves through a medium.
- **Navier-Stokes Equations**: Governs the motion of viscous fluid substances.
- **Poisson's Equation**: Used in electrostatics, mechanical engineering, and theoretical physics.

## Numerical Methods

- **Finite Difference Method (FDM)**: A straightforward approach to discretize PDEs on regular grids.
- **Finite Element Method (FEM)**: Ideal for complex geometries and adaptive mesh refinement.
- **Spectral Methods**: Provides high accuracy for smooth solutions and periodic boundary conditions.

## Usage (Already adjusted for not released v0.0.95, use Manual for accurate usage)
**Using the Package**:
```julia
    using PDEngine

    # Example: Solving the heat equation in 1D
    Δx = 0.01
    Δt = 0.0001
    T = 0.1
    α = 0.01
    N = 100
    temperature_distribution = heat(N, α, T, Δx, Δt) # solving the heat equation with the default (spectral method)

    # Example: Solving Poisson's equation in 2D
    f = zeros(N+1, N+1)  # Example source term
    Δx = 0.01
    poisson_solution = poisson_fem(N, f, Δx)  # solving the poisson equations,set to use with the finite elements method
```

## Currently working on:
### v0.0.95
- Better Function Naming
- Impementing Spectral methods
- Setting a default method with very simple Name like heat(params)
### v0.1.0 (Full release)
- Parallel computing to get the full speed of julia
- more methods and or PEDs
- ensure all Methods work without any flaws
- ...

## Papers used:

- [Lehman](https://www.lehman.edu/faculty/dgaranin/Mathematical_Physics/Mathematical_physics-13-Partial_differential_equations.pdf) : Numerical Solutions for PDEs using Mathematica
- [Spectral Methods](https://www.mech.kth.se/~ardeshir/courses/literature/Notes_Spectral_Methods.pdf) : Spectral Methods and tests by Phillip Schlatter
- [ByJus](https://byjus.com/maths/partial-differential-equation/) : Representation of a PDE
- [ScholarPedia](http://www.scholarpedia.org/article/Partial_differential_equation) : Introduction to PDEs
    
