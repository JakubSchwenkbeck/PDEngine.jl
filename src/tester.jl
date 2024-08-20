include("PDE_lib.jl")
using .PDEngine
using LinearAlgebra

function test_heat_equation(N_values, α, T, Δt)
    for N in N_values
        Δx = 1.0 / N    # Spatial step size
        
        # Analytical solution at final time T
        analytical_solution(x, t) = exp(-α * π^2 * t) * sin(π * x)

        # Solve using FDM
        u_fdm = heat_eq_1d_fdm(N, α, T, Δx, Δt)

        # Solve using FEM
        u_fem = heat_eq_1d_fem(N, α, T, Δx, Δt)

        # Compare with analytical solution
        x = 0:Δx:1
        u_exact = analytical_solution.(x, T)

        # Calculate relative error norms
        error_fdm = norm(u_fdm - u_exact) / norm(u_exact)
        error_fem = norm(u_fem - u_exact) / norm(u_exact)

        println("Grid Points: $N")
        println("Relative Error in FDM: $error_fdm")
        println("Relative Error in FEM: $error_fem")
        println("------------------------------------------------")
    end
end

# Example test parameters
N_values = [10, 20, 50, 100]  # Different grid resolutions
α = 1.0                       # Thermal diffusivity
T = 0.01                      # Total simulation time
Δt = 0.0001                   # Temporal step size

test_heat_equation(N_values, α, T, Δt)
