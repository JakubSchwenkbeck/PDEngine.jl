using LinearAlgebra,FFTW
"""
heat(N, α, T, Δx, Δt) is the default function using ...


"""

"""
    heat_spectral(N, α, T, Δx, Δt)

Solve the 1D heat equation using the Spectral Method.

# Arguments
- `N::Int`: The number of spatial grid points (excluding boundary points).
- `α::Float64`: The thermal diffusivity or the diffusion coefficient.
- `T::Float64`: The total simulation time.
- `Δx::Float64`: The spatial step size.
- `Δt::Float64`: The time step size.

# Returns
- `Vector{Float64}`: The temperature distribution at the final time step.
"""

function heat_spectral(N, α, T, Δx, Δt)
    # Define the spatial domain
    x = 0:Δx:1
    L = 1.0
    dx = L / N
    x = dx * (0:N)
    
    # Initialize the temperature field with a Gaussian profile
    u = exp.(-((x .- 0.5).^2) / (2 * (0.1)^2))
    
    # Calculate the wavenumbers (k) for the spectral method
    k = 2 * π * (0:(N ÷ 2)) / L
    
    # Precompute the Fourier transform of the initial condition
    u_hat = rfft(u)
    
    # Time-stepping loop
    num_steps = Int(T / Δt)
    for _ in 1:num_steps
        # Compute the time evolution in Fourier space
        u_hat .= u_hat .* exp.(-α * k.^2 * Δt)
    end
    
    # Inverse Fourier transform to get the solution in physical space
    u = irfft(u_hat, N + 1)
    
    return u
end





"""
    check_stability(α, Δx, Δt)

Check the stability condition for a finite difference scheme solving the heat equation.

# Arguments
- `α::Float64`: The thermal diffusivity or the diffusion coefficient.
- `Δx::Float64`: The spatial step size.
- `Δt::Float64`: The time step size.

# Returns
- `Bool`: `true` if the stability condition is satisfied, otherwise `false`.
- `String`: A message indicating whether the parameters satisfy the stability condition.
"""
function check_stability(α, Δx, Δt)
    stability_param = α * Δt / Δx^2
    if stability_param <= 0.5
        return true, "The parameters satisfy the stability condition."
    else
        return false, "The parameters do NOT satisfy the stability condition."
    end
end


"""
    heat_fdm(N, α, T, Δx, Δt)

Solve the 1D heat equation using the Finite Difference Method (FDM).

# Arguments
- `N::Int`: The number of spatial grid points (excluding boundary points).
- `α::Float64`: The thermal diffusivity or the diffusion coefficient.
- `T::Float64`: The total simulation time.
- `Δx::Float64`: The spatial step size.
- `Δt::Float64`: The time step size.

# Returns
- `Vector{Float64}`: The temperature distribution at the final time step.
"""
function heat_fdm(N, α, T, Δx, Δt)
    # Check if the stability condition is satisfied
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    # Calculate the number of time steps
    num_steps = Int(T / Δt)
    
    # Define the spatial domain
    x = 0:Δx:1
    
    # Initialize the temperature field
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set the initial condition: Gaussian profile centered at x = 0.5
    for i in 1:N+1
        u[i] = exp(-((x[i] - 0.5)^2) / (2 * (0.1)^2))
    end
    
    # Time-stepping loop
    for _ in 1:num_steps
        # Update the temperature field using the finite difference scheme
        for i in 2:N
            u_new[i] = u[i] + α * Δt / Δx^2 * (u[i+1] - 2*u[i] + u[i-1])
        end
        # Apply Dirichlet boundary conditions (u = 0 at the boundaries)
        u_new[1] = 0
        u_new[end] = 0
        
        # Swap the references to avoid unnecessary copying
        u, u_new = u_new, u
    end
    
    return u
end

"""
    heat_fem(N, α, T, Δx, Δt)

Solve the 1D heat equation using the Finite Element Method (FEM).

# Arguments
- `N::Int`: The number of spatial grid points (excluding boundary points).
- `α::Float64`: The thermal diffusivity or the diffusion coefficient.
- `T::Float64`: The total simulation time.
- `Δx::Float64`: The spatial step size.
- `Δt::Float64`: The time step size.

# Returns
- `Vector{Float64}`: The temperature distribution at the final time step.
"""
function heat_fem(N, α, T, Δx, Δt)
    # Check if the stability condition is satisfied
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    # Calculate the number of time steps
    num_steps = Int(T / Δt)
    
    # Define the spatial domain (nodes)
    nodes = 0:Δx:1
    
    # Initialize the temperature field
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set the initial condition: heat source at the middle of the domain
    mid_point = round(Int, N/2) + 1
    u[mid_point] = 1.0
    
    # Initialize the stiffness (K) and mass (M) matrices
    K = zeros(N+1, N+1)
    M = zeros(N+1, N+1)
    
    # Assemble the stiffness and mass matrices
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
    
    # Apply Dirichlet boundary conditions (u = 0 at the boundaries)
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
        
        # Apply boundary conditions to the load vector
        f[1] = 0
        f[end] = 0
        
        # Solve for the new temperature field
        u_new = (M + α * Δt / 2 * K) \ (M * u + α * Δt / 2 * f)
        
        # Update the temperature field for the next time step
        u[:] = u_new[:]
    end
    
    return u
end

"""
    heat_crank_nicolson(N, α, T, Δx, Δt)

Solve the 1D heat equation using the Crank-Nicolson method.

# Arguments
- `N::Int`: The number of spatial grid points (excluding boundary points).
- `α::Float64`: The thermal diffusivity or the diffusion coefficient.
- `T::Float64`: The total simulation time.
- `Δx::Float64`: The spatial step size.
- `Δt::Float64`: The time step size.

# Returns
- `Vector{Float64}`: The temperature distribution at the final time step.
"""
function heat_nicolson(N, α, T, Δx, Δt)
    # Check if the stability condition is satisfied
    stability, message = check_stability(α, Δx, Δt)
    if !stability
        println(message)
    end
    
    # Calculate the number of time steps
    num_steps = Int(T / Δt)
    
    # Define the spatial domain
    x = 0:Δx:1
    
    # Initialize the temperature field
    u = zeros(N+1)
    u_new = similar(u)
    
    # Set the initial condition: heat source at the middle of the domain
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
    
    # Apply Dirichlet boundary conditions (u = 0 at the boundaries)
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
        
        # Update the temperature field for the next time step
        u[:] = u_new[:]
    end
    
    return u
end
