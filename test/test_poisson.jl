

include("../src/PDE_lib.jl")
using .PDEngine
using LinearAlgebra


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