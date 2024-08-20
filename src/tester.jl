# Testing code
include("PDE_lib.jl")
using .PDEngine

using LinearAlgebra

function test_against_analytical()
    N = 2          # Number of spatial grid points
    α = 1.0         # Thermal diffusivity
    T = 0.01        # Total simulation time
    Δx = 1.0 / N    # Spatial step size
    Δt = 0.0001     # Temporal step size

    # Analytical solution at final time T
    analytical_solution(x, t) = exp(-α * π^2 * t) * sin(π * x)

    # Solve using FDM
    u_fdm = heat_eq_1d_fdm(N, α, T, Δx, Δt)

    # Solve using FEM
    u_fem = heat_eq_1d_fem(N, α, T, Δx, Δt)

    # Compare with analytical solution
    x = 0:Δx:1
    u_exact = analytical_solution.(x, T)

    # Calculate error norms
    error_fdm = norm(u_fdm - u_exact)
    error_fem = norm(u_fem - u_exact)

    println("Error in FDM: $error_fdm")
    println("Error in FEM: $error_fem")
end

test_against_analytical()
    