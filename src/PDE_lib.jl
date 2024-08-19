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

"""
Solve the 2D Navier-Stokes equations using finite difference method.
The Navier-Stokes equations describe the motion of viscous fluid substances.
"""
function navier_stokes_2d(N, ν, T, Δx, Δt)
    # Initialize velocity and pressure fields
    u = zeros(N+1, N+1)   # x-component of velocity
    v = zeros(N+1, N+1)   # y-component of velocity
    p = zeros(N+1, N+1)   # Pressure field
    u_new, v_new = similar(u), similar(v)

    # Time-stepping loop
    for t in Δt:Δt:T
        # Update velocity and pressure fields using finite difference methods
        # This section is a placeholder; actual implementation requires solving complex equations
        
        # Apply boundary conditions and pressure correction
        # (Implement pressure-correction method like SIMPLE or PISO)
    end

    return u, v, p
end

"""
Solve 2D Poisson's equation using the finite difference method.
Poisson's equation is given by: ∇²u = -f
where f is a source term.
"""
function poisson_2d(N, f, Δx)
    # Initialize solution field
    u = zeros(N+1, N+1)
    u_new = similar(u)

    # Iterative solver to update the solution field
    for _ in 1:1000  # Number of iterations (could be adjusted or replaced with convergence check)
        # Update the solution field using finite difference scheme
        for i in 2:N
            for j in 2:N
                u_new[i,j] = 0.25 * (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - Δx^2 * f[i,j])
            end
        end
        
        # Swap references to avoid unnecessary copying
        u[:] = u_new[:]
    end

    return u
end
