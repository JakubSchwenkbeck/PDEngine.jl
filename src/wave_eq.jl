

"""
Solve the 1D wave equation using the finite difference method.
The wave equation is given by: ∂²u/∂t² = c² ∇²u
where c is the wave speed.
"""
function wave_fdm(N, c, T, Δx, Δt)
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
Solve the 1D wave equation using the finite element method (FEM).
The wave equation is given by: ∂²u/∂t² = c² ∇²u
where c is the wave speed.

Arguments:
- N: Number of spatial grid points (excluding boundary points)
- c: Wave speed
- T: Total simulation time
- Δx: Spatial step size
- Δt: Temporal step size

Returns:
- u: The final displacement field
"""
function wave_fem(N, c, T, Δx, Δt)
    # Number of time steps
    num_steps = Int(T / Δt)
    
    # Create a grid of spatial points
    nodes = 0:Δx:1

    # Initialize displacement and velocity fields
    u = zeros(N+1)     # Displacement
    u_prev = similar(u) # Displacement at previous time step
    u_next = similar(u) # Displacement at next time step

    # Initialize stiffness matrix (K) and mass matrix (M)
    K = zeros(N+1, N+1)
    M = zeros(N+1, N+1)

    # Assemble stiffness matrix and mass matrix
    for i in 1:N
        # Local matrices for element i
        x1 = nodes[i]
        x2 = nodes[i+1]
        L = x2 - x1
        k_local = (c^2 / L) * [1 -1; -1 1]
        m_local = L / 6 * [2 1; 1 2]

        # Assemble into global stiffness and mass matrices
        K[i:i+1, i:i+1] += k_local
        M[i:i+1, i:i+1] += m_local
    end

    # Apply boundary conditions (Dirichlet: u = 0 at the boundaries)
    K[1, :] .= 0
    K[1, 1] = 1
    M[1, :] .= 0
    M[1, 1] = 1
    K[end, :] .= 0
    K[end, end] = 1
    M[end, :] .= 0
    M[end, end] = 1

    # Time-stepping loop
    for _ in 1:num_steps
        # Compute the right-hand side vector
        f = M * u_next - K * u_prev
        
        # Solve for the new displacement field
        u_next = K \ f
        
        # Update displacement fields
        u_prev, u, u_next = u, u_next, u_prev
    end

    return u
end
