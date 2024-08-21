# PDE Library Manual

Welcome to the manual for the PDEngine Library. This document provides detailed instructions on how to use the functions available in the library for solving various partial differential equations (PDEs) using Julia.

## Heat Equation Solver

### Finite Difference Method (FDM)

**Function Signature:**  
`function heat_eq_1d_fdm(N, α, T, Δx, Δt)`

**Description:**  
Solves the 1D heat equation using the finite difference method.

**Arguments:**
- `N`: Number of grid points (excluding boundaries)
- `α`: Thermal diffusivity
- `T`: Total simulation time
- `Δx`: Spatial step size
- `Δt`: Temporal step size

**Returns:**
- `u`: Final temperature distribution

**Example Usage:**  
To solve the 1D heat equation using FDM with a grid of 50 points, thermal diffusivity of 1, and a total simulation time of 1.0, use the following parameters:
- `N = 50`
- `α = 1`
- `T = 1.0`
- `Δx = 1.0 / N`
- `Δt = 0.001`

Call the function as follows:
```julia
u_fdm = heat_eq_1d_fdm(N, α, T, Δx, Δt)
```

Finite Element Method (FEM)
Function Signature:
function heat_eq_1d_fem(N, α, T, Δx, Δt)

Description:
Solves the 1D heat equation using the finite element method.

Arguments:

N: Number of grid points (excluding boundaries)
α: Thermal diffusivity
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: Final temperature distribution
Example Usage:
To solve the 1D heat equation using FEM with the same parameters as above, use:

N = 50
α = 1
T = 1.0
Δx = 1.0 / N
Δt = 0.001
Call the function as follows:

```julia

u_fem = heat_eq_1d_fem(N, α, T, Δx, Δt)
```

Crank-Nicolson Method (CN)
Function Signature:
function heat_eq_1d_crank_nicolson(N, α, T, Δx, Δt)

Description:
Solves the 1D heat equation using the Crank-Nicolson method.

Arguments:

N: Number of grid points (excluding boundaries)
α: Thermal diffusivity
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: Final temperature distribution
Example Usage:
To solve the 1D heat equation using CN with the same parameters, use:

N = 50
α = 1
T = 1.0
Δx = 1.0 / N
Δt = 0.001
Call the function as follows:

```julia

u_cn = heat_eq_1d_crank_nicolson(N, α, T, Δx, Δt)
```
Navier-Stokes Solver
2D Navier-Stokes Equations
Function Signature:
function navier_stokes_2d_solver(N, Re, T, Δx, Δt)

Description:
Solves the 2D Navier-Stokes equations for incompressible flow.

Arguments:

N: Number of grid points in each dimension (excluding boundaries)
Re: Reynolds number
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: Velocity field
p: Pressure field
Example Usage:
To solve the 2D Navier-Stokes equations with a grid of 50x50 points and a Reynolds number of 100, use:

N = 50
Re = 100
T = 1.0
Δx = 1.0 / N
Δt = 0.001
Call the function as follows:

```julia

u, p = navier_stokes_2d_solver(N, Re, T, Δx, Δt)
```
Poisson's Equation Solver
Finite Difference Method (FDM)
Function Signature:
function poisson_2d_fdm(N, f, Δx)

Description:
Solves 2D Poisson's equation using the finite difference method. Poisson's equation is given by: ∇²u = -f where f is a source term.

Arguments:

N: Number of grid points along each dimension (excluding boundaries)
f: Source term (function or matrix)
Δx: Spatial step size
Returns:

u: The final solution field
Example Usage:
To solve 2D Poisson's equation with a source term and a grid of 50x50 points:

N = 50
f = ... (define your source term)
Δx = 1.0 / N
Call the function as follows:

```julia

u = poisson_2d_fdm(N, f, Δx)
```
Finite Element Method (FEM)
Function Signature:
function poisson_2d_fem(N, f, Δx)

Description:
Solves 2D Poisson's equation using the finite element method (FEM). Poisson's equation is given by: ∇²u = -f where f is a source term.

Arguments:

N: Number of grid points along each dimension (excluding boundaries)
f: Source term (function)
Δx: Spatial step size
Returns:

u: The final solution field in grid format
Example Usage:
To solve 2D Poisson's equation using FEM with a grid of 50x50 points:

N = 50
f = ... (define your source term)
Δx = 1.0 / N
Call the function as follows:

```julia

u = poisson_2d_fem(N, f, Δx)
```
Wave Equation Solver
Finite Difference Method (FDM)
Function Signature:
function wave_eq_1d_fdm(N, c, T, Δx, Δt)

Description:
Solves the 1D wave equation using the finite difference method. The wave equation is given by: ∂²u/∂t² = c² ∇²u where c is the wave speed.

Arguments:

N: Number of spatial grid points (excluding boundary points)
c: Wave speed
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: The final displacement field
Example Usage:
To solve the 1D wave equation with a grid of 50 points and a wave speed of 1.0:

N = 50
c = 1.0
T = 1.0
Δx = 1.0 / N
Δt = 0.001
Call the function as follows:

```julia

u = wave_eq_1d_fdm(N, c, T, Δx, Δt)
```
Finite Element Method (FEM)
Function Signature:
function wave_eq_1d_fem(N, c, T, Δx, Δt)

Description:
Solves the 1D wave equation using the finite element method (FEM). The wave equation is given by: ∂²u/∂t² = c² ∇²u where c is the wave speed.

Arguments:

N: Number of spatial grid points (excluding boundary points)
c: Wave speed
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: The final displacement field
Example Usage:
To solve the 1D wave equation using FEM with a grid of 50 points:

N = 50
c = 1.0
T = 1.0
Δx = 1.0 / N
Δt = 0.001
Call the function as follows:

```julia

u = wave_eq_1d_fem(N, c, T, Δx, Δt)```