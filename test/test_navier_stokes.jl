using Test
using LinearAlgebra
using Statistics


# Import the navier_stokes_2d function from the source file
include("../src/navier_stokes_eq.jl")

# Test Suite for the 2D Navier-Stokes Solver

# Test 1: Check if the function runs without errors and returns correct output types
@testset "Navier-Stokes 2D Solver Tests" begin
    # Test parameters
    N = 50
    ν = 0.01
    T = 1.0
    Δx = 1.0 / N
    Δt = 0.001

    # Run the solver
    u, v, p = navier_stokes_2d(N, ν, T, Δx, Δt)

    # Test 1.1: Ensure the output is the correct size
    @test size(u) == (N+1, N+1)
    @test size(v) == (N+1, N+1)
    @test size(p) == (N+1, N+1)

    # Test 1.2: Ensure the output values are finite numbers
    @test all(isfinite, u)
    @test all(isfinite, v)
    @test all(isfinite, p)

    # Test 2: Boundary conditions - check if boundary values are correctly set
    @test all(u[1, :] .== 0)       # u = 0 at y = 0
    @test all(u[end, :] .== 0)     # u = 0 at y = 1
    @test all(v[:, 1] .== 0)       # v = 0 at x = 0
    @test all(v[:, end] .== 0)     # v = 0 at x = 1

    # Test 3: Check pressure field boundary conditions
    @test isapprox(mean(p[:, 1]), mean(p[:, 2]), atol=1e-6)
    @test isapprox(mean(p[:, end]), mean(p[:, end-1]), atol=1e-6)
    @test isapprox(mean(p[1, :]), mean(p[2, :]), atol=1e-6)
    @test isapprox(mean(p[end, :]), mean(p[end-1, :]), atol=1e-6)

    # Test 4: Consistency Check
    # We expect the velocity fields to be very small due to the small timestep and viscosity
    @test maximum(abs.(u)) < 1e-2
    @test maximum(abs.(v)) < 1e-2
end


