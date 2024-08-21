# PDE Library Reference

This document provides detailed API documentation for the functions in the PDE Library. Each function is described with its signature, arguments, and return values.

## Functions

### `heat_eq_1d_fdm`

**Signature:**

```julia
function heat_eq_1d_fdm(N, α, T, Δx, Δt)
```
Description:

Solves the 1D heat equation using the finite difference method.

Arguments:

N: Number of grid points (excluding boundaries)
α: Thermal diffusivity
T: Total simulation time
Δx: Spatial step size
Δt: Temporal step size
Returns:

u: Final temperature distribution
heat_eq_1d_fem
Signature:

```julia

function heat_eq_1d_fem(N, α, T, Δx, Δt)
```
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
heat_eq_1d_crank_nicolson
Signature:

```julia

function heat_eq_1d_crank_nicolson(N, α, T, Δx, Δt)
```
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
navier_stokes_2d_solver
Signature:

```julia

function navier_stokes_2d_solver(N, Re, T, Δx, Δt)
```
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
poisson_2d_fdm
Signature:

```julia

function poisson_2d_fdm(N, f, Δx)
```
Description:

Solves 2D Poisson's equation using the finite difference method. Poisson's equation is given by: ∇²u = -f where f is a source term.

Arguments:

N: Number of grid points along each dimension (excluding boundaries)
f: Source term (function or matrix)
Δx: Spatial step size
Returns:

u: The final solution field
poisson_2d_fem
Signature:

```julia

function poisson_2d_fem(N, f, Δx)
```
Description:

Solves 2D Poisson's equation using the finite element method (FEM). Poisson's equation is given by: ∇²u = -f where f is a source term.

Arguments:

N: Number of grid points along each dimension (excluding boundaries)
f: Source term (function)
Δx: Spatial step size
Returns:

u: The final solution field in grid format
wave_eq_1d_fdm
Signature:

```julia

function wave_eq_1d_fdm(N, c, T, Δx, Δt)
```
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
wave_eq_1d_fem
Signature:

```julia

function wave_eq_1d_fem(N, c, T, Δx, Δt)
```
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
