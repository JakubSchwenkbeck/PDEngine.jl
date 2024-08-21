using Test
using LinearAlgebra
include("../src/PDE_lib.jl")  
using .PDEngine

function test_heat_equation()
    # Test parameters
    heat_N_values = [5, 10, 20, 50, 100, 300, 500]  # Grid resolutions
    heat_α = 1                        
    heat_T = 1e-3                      # Total time
    heat_Δt = 1e-4                     # Time step size

    @testset "Heat Equation Solver Tests" begin
        for N in heat_N_values
            Δx = 1.0 / N    # Spatial step size
            
            # Analytical solution at final time T
            analytical_solution(x, t) = exp(-heat_α * π^2 * t) * sin(π * x)

            # Grid points including both boundaries (0 and 1)
            x = 0:Δx:1
            
            # Solve using FDM
            u_fdm = heat_eq_1d_fdm(N, heat_α, heat_T, Δx, heat_Δt)

            # Solve using FEM
            u_fem = heat_eq_1d_fem(N, heat_α, heat_T, Δx, heat_Δt)

            # Solve using Crank-Nicolson
            u_cn = heat_eq_1d_crank_nicolson(N, heat_α, heat_T, Δx, heat_Δt)

            # Analytical solution at final time
            u_exact = analytical_solution.(x, heat_T)

            # Ensure that the numerical and analytical solutions have the same length
            @test length(u_fdm) == length(u_exact) "FDM result length mismatch"
            @test length(u_fem) == length(u_exact) "FEM result length mismatch"
            @test length(u_cn) == length(u_exact) "CN result length mismatch"

            # Calculate relative error norms
            error_fdm = abs((norm(u_fdm - u_exact) / norm(u_exact)) - 1)
            error_fem = abs((norm(u_fem - u_exact) / norm(u_exact)) - 1)
            error_cn = abs((norm(u_cn - u_exact) / norm(u_exact)) - 1)

            # Relative error checks
            @test error_fdm < 1e-2 "High relative error in FDM for N=$N"
            @test error_fem < 1e-2 "High relative error in FEM for N=$N"
            @test error_cn < 1e-2 "High relative error in CN for N=$N"
        end
    end
end

# Run the test
test_heat_equation()
