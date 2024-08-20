include("PDE_lib.jl")
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

#test_heat_equation(heat_N_values, heat_α, heat_T, heat_Δt)










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

#test_wave_equation(wave_N_values, wave_c, wave_T, wave_Δt)


"""
Test the 2D Poisson equation solver 
Compares the numerical solution against an analytical solution.
"""


function test_poisson_2d(N_values)
    # Parameters
    for N in N_values
    Δx = 1.0 / N
    
    # Define the source term f as a constant function for simplicity
    f(x, y) = 1.0

    # Convert to matrix form for finite difference method
    f_matrix = [f(x, y) for x in 0:Δx:1, y in 0:Δx:1]

    # Solve using FDM
    u_fdm = poisson_2d_fdm(N, f_matrix, Δx)

    # Solve using FEM
    f_func(x, y) = 1.0
    u_fem = poisson_2d_fem(N, f_func, Δx)

    # Calculate the exact solution (for comparison, assuming f=1 leads to a simple linear solution)
    u_exact = [x * (1 - x) * y * (1 - y) for x in 0:Δx:1, y in 0:Δx:1]

    # Calculate the relative error norms
    error_fdm = norm(u_fdm - u_exact) / norm(u_exact)
    error_fem = norm(u_fem - u_exact) / norm(u_exact)

    println("Grid Points: $N")
    println("Relative Error in FDM: $error_fdm")
    println("Relative Error in FEM: $error_fem")
    println("------------------------------------------------")
    end
end

poisson_N_values = [5, 10, 20, 50, 100]
test_poisson_2d(poisson_N_values)