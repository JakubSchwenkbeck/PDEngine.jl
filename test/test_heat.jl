include("../src/PDE_lib.jl")
using .PDEngine
using LinearAlgebra

function test_heat_equation(N_values, α, T, Δt)
    println("Heat Equation Tester : \n")
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
heat_N_values = [5,10, 20, 50, 100,500,1000]  # Different grid resolutions
heat_α = 1.0                       # Thermal diffusivity
heat_T = 1e-10                    # Total simulation time --v
heat_Δt = 1e-10            # Temporal step size      --> Biggest impact on accuracy!

test_heat_equation(heat_N_values, heat_α, heat_T, heat_Δt)
