
using LinearAlgebra
include("../src/PDEngine.jl")  
using .PDEngine

function test_heat_equation()
    # Test parameters
    heat_N_values = [5, 10, 20, 50, 100, 300, 500]  # Grid resolutions
    heat_α = 1                        
    heat_T = 1e-3                      # Total time
    heat_Δt = 1e-4                     # Time step size

 
        for N in heat_N_values
            Δx = 1.0 / N    # Spatial step size
            
            # Analytical solution at final time T
            analytical_solution(x, t) = exp(-heat_α * π^2 * t) * sin(π * x)

            # Grid points including both boundaries (0 and 1)
            x = 0:Δx:1
            
            # Solve using FDM
            u_fdm = heat_fdm(N, heat_α, heat_T, Δx, heat_Δt)

            # Solve using FEM
            u_fem = heat_fem(N, heat_α, heat_T, Δx, heat_Δt)

            # Solve using Crank-Nicolson
            u_cn = heat_nicolson(N, heat_α, heat_T, Δx, heat_Δt)

            u_sp = heat_spectral(N, heat_α, heat_T, heat_Δt)
            # Analytical solution at final time
            u_exact = analytical_solution.(x, heat_T)

            # Ensure that the numerical and analytical solutions have the same length
          
            # Calculate relative error norms
            error_fdm = abs((norm(u_fdm - u_exact) / norm(u_exact)) - 0)
            error_fem = abs((norm(u_fem - u_exact) / norm(u_exact)) - 0)
            error_cn = abs((norm(u_cn - u_exact) / norm(u_exact)) - 0)
            error_sp = abs((norm(u_sp - u_exact) / norm(u_exact)) - 0)

            # Relative error checks
            if error_fdm < 1e-2
                 print("High relative error in FDM for N=$N" )
            end
            if error_fem < 1e-2 
               print( "High relative error in FEM for N=$N")
            end
            if error_cn < 1e-2 
                print("High relative error in CN for N=$N")
            end
            if error_sp < 1e-2 
                print("High relative error in SP for N=$N")
            end
            println("Grid Points: $N")
            println("Relative Error in FDM: $error_fdm")
            println("Relative Error in FEM: $error_fem")
            println("Relative Error in CN: $error_cn")
            println("Relative Error in SP: $error_sp")
            println("------------------------------------------------")
        end
    end

# Run the test
test_heat_equation()
