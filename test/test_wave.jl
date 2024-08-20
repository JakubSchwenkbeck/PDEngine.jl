
include("../src/PDE_lib.jl")
using .PDEngine
using LinearAlgebra

function test_wave_equation(N_values, c, T, Δt)
    println("Wave Equation Tester : \n")
    for N in N_values
        Δx = 1.0 / N    # Spatial step size

        # Analytical solution at final time T
        analytical_solution(x, t) = sin(π * x) * cos(c * π * t)

        # Solve using FDM
        u_fdm = wave_eq_1d_fdm(N, c, T, Δx, Δt)

        # Solve using FEM
        u_fem = wave_eq_1d_fem(N, c, T, Δx, Δt)

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
wave_N_values = [5, 10, 20, 50, 100, 500, 1000]  # Different grid resolutions
wave_c = 1e-4                        # Wave speed
wave_T = 1e-7                        # Total simulation time (ensure it is an appropriate value for the wave propagation)
wave_Δt = 1e-7              # Temporal step size

test_wave_equation(wave_N_values, wave_c, wave_T, wave_Δt)
