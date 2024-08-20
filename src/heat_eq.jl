"""
Solve the 1D heat equation using the finite difference method (FDM).
The heat equation is given by: ∂u/∂t = α ∇²u
where α is the thermal diffusivity.
"""
function heat_eq_1d_fdm(N, α, T, Δx, Δt)
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
          # Apply boundary conditions (Dirichlet: u = 0 at the boundaries)
          u_new[1] = 0
          u_new[end] = 0
                 
        # Swap references to avoid unnecessary copying
        u, u_new = u_new, u
    end

    return u
end


"""
Solve the 1D heat equation using the finite element method (FEM).
The heat equation is given by: ∂u/∂t = α ∇²u
where α is the thermal diffusivity.

Arguments:
- N: Number of spatial grid points (excluding boundary points)
- α: Thermal diffusivity
- T: Total simulation time
- Δx: Spatial step size
- Δt: Temporal step size

Returns:
- u: The final temperature field
"""
function heat_eq_1d_fem(N, α, T, Δx, Δt)
    # Number of time steps
    num_steps = Int(T / Δt)
    
    # Create a grid of spatial points
    nodes = 0:Δx:1

    # Initialize temperature field and new temperature field
    u = zeros(N+1)
    u_new = similar(u)

    # Set initial condition: a heat source in the middle of the domain
    mid_point = round(Int, N/2) + 1  # Adjusted for 1-based indexing in Julia
    u[mid_point] = 1.0

    # Initialize stiffness matrix (K) and load vector (f)
    K = zeros(N+1, N+1)
    f = zeros(N+1)

    # Assemble stiffness matrix and load vector
    for i in 1:N
        # Local stiffness matrix for element i
        x1 = nodes[i]
        x2 = nodes[i+1]
        L = x2 - x1
        k_local = α / L * [1 -1; -1 1]

        # Assemble into global stiffness matrix
        K[i:i+1, i:i+1] += k_local
    end

    # Apply boundary conditions (Dirichlet: u = 0 at the boundaries)
    K[1, :] .= 0
    K[1, 1] = 1
    f[1] = 0
    K[end, :] .= 0
    K[end, end] = 1
    f[end] = 0

    # Time-stepping loop
    for _ in 1:num_steps
        # Update the load vector based on the previous solution
        f[2:end-1] = u[2:end-1]  # HERE, no external source term!!!
        
        # Solve for the new temperature field
        u_new = K \ f
        
        # Update temperature field for the next time step
        u[:] = u_new[:]
    end

    return u
end
