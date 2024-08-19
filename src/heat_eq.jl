"""
Solve the 1D heat equation using the finite difference method.
The heat equation is given by: ∂u/∂t = α ∇²u
where α is the thermal diffusivity.
"""
function heat_eq_1d(N, α, T, Δx, Δt)
    # Create a grid of spatial points
    x = 0:Δx:1

    # Initialize temperature field and new temperature field
    u = zeros(N+1)
    u_new = similar(u)

    # Set initial condition: a heat source in the middle of the domain
    u[round(Int, N/2)] = 1.0

    # Time-stepping loop
    for t in 0:Δt:T
        # Update the temperature field using finite difference scheme
        for i in 2:N
            u_new[i] = u[i] + α * Δt / Δx^2 * (u[i+1] - 2*u[i] + u[i-1])
        end
        
        # Swap references to avoid unnecessary copying
        u, u_new = u_new, u
    end

    return u
end
