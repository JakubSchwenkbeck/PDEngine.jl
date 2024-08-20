include("../src/PDE_lib.jl")
using .PDEngine
using LinearAlgebra

function test_heat_equation(N_values, α, T, Δt)
    println("Heat Equation Tester : \n")
    for N in N_values
        Δx = 1.0 / N    # Spatial step size
        
        # Analytical solution at final time T
        analytical_solution(x, t) = exp(-α * π^2 * t) * sin(π * x)

        # Grid points including both boundaries (0 and 1)
        x = 0:Δx:1
        
        # Solve using FDM
        u_fdm = heat_eq_1d_fdm(N, α, T, Δx, Δt)

        # Solve using FEM
        u_fem = heat_eq_1d_fem(N, α, T, Δx, Δt)

        u_cn = heat_eq_1d_crank_nicolson(N, α, T, Δx, Δt)
        # Compare with analytical solution
        u_exact = analytical_solution.(x, T)

        # Ensure that u_fdm, u_fem, and u_exact have the same length
        if length(u_fdm) != length(u_exact)
            println("Error: FDM result length mismatch with analytical solution.")
        end
        if length(u_fem) != length(u_exact)
            println("Error: FEM result length mismatch with analytical solution.")
        end

        # Calculate relative error norms
        error_fdm = abs((norm(u_fdm - u_exact) / norm(u_exact) )-1)
        error_fem = abs(norm((u_fem - u_exact) / norm(u_exact))-1)
        error_cn = abs((norm(u_cn - u_exact) / norm(u_exact))-1)
        println("Grid Points: $N")
        println("Relative Error in FDM: $error_fdm")
        println("Relative Error in FEM: $error_fem")
        println("Relative Error in CN: $error_cn")
        println("-------------------------------------------------")
    end
end

# Example test parameters
heat_N_values = [5, 10, 20, 50, 100,300,500]  # Reduced grid resolutions for testing
heat_α = 1                         
heat_T = 1e-3                        # Increased time for better evaluation
heat_Δt = 1e-4                        # Smaller time step for accuracy

test_heat_equation(heat_N_values, heat_α, heat_T, heat_Δt)