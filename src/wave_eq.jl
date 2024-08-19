

"""
Solve the 1D wave equation using the finite difference method.
The wave equation is given by: ∂²u/∂t² = c² ∇²u
where c is the wave speed.
"""
function wave_eq_1d(N, c, T, Δx, Δt)
    # Create a grid of spatial points
    x = 0:Δx:1

    # Initialize velocity and displacement fields
    u = zeros(N+1)     # Displacement
    v = zeros(N+1)     # Velocity (time derivative of displacement)
    u_prev = similar(u)
    u_next = similar(u)

    # Set initial condition: a displacement in the middle of the domain
    u[round(Int, N/2)] = 1.0
    u_prev = u

    # Time-stepping loop
    for t in Δt:Δt:T
        # Update the displacement field using finite difference scheme
        for i in 2:N
            u_next[i] = 2 * u[i] - u_prev[i] + (c^2 * Δt^2 / Δx^2) * (u[i+1] - 2 * u[i] + u[i-1])
        end
        
        # Swap references to avoid unnecessary copying
        u_prev, u, u_next = u, u_next, u_prev
    end

    return u
end