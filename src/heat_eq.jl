using LinearAlgebra

# Stability Check Function
function check_stability(α, Δx, Δt)
    stability_param = α * Δt / Δx^2
    if stability_param <= 0.5
        return true, "The parameters satisfy the stability condition."
    else
        return false, "The parameters do NOT satisfy the stability condition."
    end
end

# Finite Difference Method (FDM)
function heat_eq_1d_fdm(N, α, T, Δx, Δt)
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    num_steps = Int(T / Δt)
    x = 0:Δx:1
    
    # Initialize temperature field
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set initial condition: Gaussian profile
    for i in 1:N+1
        u[i] = exp(-((x[i] - 0.5)^2) / (2 * (0.1)^2))
    end
    
    # Time-stepping loop
    for _ in 1:num_steps
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

# Finite Element Method (FEM)
function heat_eq_1d_fem(N, α, T, Δx, Δt)
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    num_steps = Int(T / Δt)
    nodes = 0:Δx:1
    
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set initial condition: a heat source in the middle of the domain
    mid_point = round(Int, N/2) + 1
    u[mid_point] = 1.0
    
    # Initialize stiffness (K) and mass (M) matrices
    K = zeros(N+1, N+1)
    M = zeros(N+1, N+1)
    
    # Assemble stiffness and mass matrices
    for i in 1:N
        x1 = nodes[i]
        x2 = nodes[i+1]
        L = x2 - x1
        k_local = α / L * [1 -1; -1 1]
        m_local = L / 6 * [2 1; 1 2]
        
        # Assemble into global matrices
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
        # Create the load vector
        f = M * u
        
        # Apply boundary conditions to load vector
        f[1] = 0
        f[end] = 0
        
        # Solve for the new temperature field
        u_new = (M + α * Δt / 2 * K) \ (M * u + α * Δt / 2 * f)
        
        # Update temperature field for the next time step
        u[:] = u_new[:]
    end
    
    return u
end

# Crank-Nicolson Method (CN)
function heat_eq_1d_crank_nicolson(N, α, T, Δx, Δt)
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    num_steps = Int(T / Δt)
    x = 0:Δx:1
    
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set initial condition: a heat source in the middle of the domain
    mid_point = round(Int, N/2) + 1
    u[mid_point] = 1.0
    
    # Coefficients for the Crank-Nicolson method
    r = α * Δt / (2 * Δx^2)
    
    # Construct the tridiagonal matrices A and B
    A = Matrix{Float64}(I, N+1, N+1)
    B = Matrix{Float64}(I, N+1, N+1)
    
    for i in 2:N
        A[i, i-1] = -r
        A[i, i] = 1 + 2 * r
        A[i, i+1] = -r
        
        B[i, i-1] = r
        B[i, i] = 1 - 2 * r
        B[i, i+1] = r
    end
    
    # Apply boundary conditions
    A[1, :] .= 0
    A[1, 1] = 1
    A[end, :] .= 0
    A[end, end] = 1
    
    B[1, :] .= 0
    B[1, 1] = 1
    B[end, :] .= 0
    B[end, end] = 1
    
    # Time-stepping loop
    for _ in 1:num_steps
        # Compute the right-hand side vector
        f = B * u
        
        # Solve the system of linear equations
        u_new = A \ f
        
        # Update temperature field for the next time step
        u[:] = u_new[:]
    end
    
    return u
end

