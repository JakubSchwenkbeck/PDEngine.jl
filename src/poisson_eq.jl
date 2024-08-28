using FFTW  # For FFT operations



function poisson(N,f,Δx)
poisson_spectral(N,f,Δx)


end
"""
Solve 2D Poisson's equation using the finite difference method.

# Arguments
- `N::Int`: Number of grid points along each dimension (excluding boundaries).
- `f::Array{Float64,2}`: A 2D array representing the source term on the grid.
- `Δx::Float64`: Spatial step size.

# Returns
- `u::Array{Float64,2}`: The final solution field.
"""
function poisson_fdm(N, f, Δx)
    # Initialize solution field
    u = zeros(N+1, N+1)
    u_new = similar(u)

    # Iterative solver to update the solution field
    for _ in 1:1000  # Number of iterations (could be adjusted or replaced with convergence check)
        # Update the solution field using finite difference scheme
        for i in 2:N
            for j in 2:N
                u_new[i, j] = 0.25 * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - Δx^2 * f[i, j])
            end
        end
        
        # Swap references to avoid unnecessary copying
        u, u_new = u_new, u
    end

    return u
end

"""
Solve 2D Poisson's equation using the finite element method (FEM).

# Arguments
- `N::Int`: Number of grid points along each dimension (excluding boundaries).
- `f::Function`: A function representing the source term.
- `Δx::Float64`: Spatial step size.

# Returns
- `u::Array{Float64,2}`: The final solution field.
"""
function poisson_fem(N, f, Δx)
    # Number of nodes
    num_nodes = (N + 1) * (N + 1)

    # Initialize matrices
    K = zeros(num_nodes, num_nodes)
    F = zeros(num_nodes)

    # Node index function
    function node_index(i, j)
        return (i - 1) * (N + 1) + j
    end

    # Loop over each element to assemble the global stiffness matrix and load vector
    for i in 1:N
        for j in 1:N
            # Local node indices (in a grid)
            n1 = node_index(i, j)
            n2 = node_index(i+1, j)
            n3 = node_index(i, j+1)
            n4 = node_index(i+1, j+1)

            # Element stiffness matrix and load vector (2x2 element)
            k_local = (1 / (4 * Δx^2)) * [2 -1 -1 0; -1 2 0 -1; -1 0 2 -1; 0 -1 -1 2]
            f_local = (Δx^2 / 4) * [f(i, j); f(i+1, j); f(i, j+1); f(i+1, j+1)]

            # Assemble into global matrices
            K[[n1, n2, n3, n4], [n1, n2, n3, n4]] .+= k_local
            F[[n1, n2, n3, n4]] .+= f_local
        end
    end

    # Apply boundary conditions (Dirichlet: u = 0 on the boundaries)
    boundary_nodes = Set()
    for i in 1:(N+1)
        push!(boundary_nodes, node_index(1, i))
        push!(boundary_nodes, node_index(N+1, i))
        push!(boundary_nodes, node_index(i, 1))
        push!(boundary_nodes, node_index(i, N+1))
    end

    for node in boundary_nodes
        K[node, :] .= 0
        K[node, node] = 1
        F[node] = 0
    end

    # Solve the linear system
    u = K \ F

    # Reshape the solution to a grid format
    u_grid = reshape(u, (N + 1, N + 1))

    return u_grid
end

"""
Solve 2D Poisson's equation using the spectral method with FFTW.

# Arguments
- `N::Int`: Number of grid points along each dimension (excluding boundaries).
- `f::Function`: A function representing the source term.
- `Δx::Float64`: Spatial step size.

# Returns
- `u::Array{Float64,2}`: The final solution field.
"""
function poisson_spectral(N, f, Δx)
    # Create grid
    x = Δx * (0:N)
    y = Δx * (0:N)
    X, Y = meshgrid(x, y)
    
    # Create source term array
    F = [f(xi, yi) for xi in x, yi in y]
    
    # Perform FFT
    F_hat = fft2(F)
    
    # Frequency arrays
    kx = 2π * [0:N÷2; -N÷2+1:-1] / (N*Δx)
    ky = 2π * [0:N÷2; -N÷2+1:-1] / (N*Δx)
    
    # Create frequency grid
    KX, KY = meshgrid(kx, ky)
    
    # Compute Laplacian in frequency domain
    L = - (KX.^2 .+ KY.^2)
    
    # Avoid division by zero (at zero frequency)
    L[L .== 0] .= 1
    
    # Solve in frequency domain
    U_hat = F_hat ./ L
    
    # Perform inverse FFT to get solution in spatial domain
    u = ifft2(U_hat)
    
    # Return real part of the solution
    return real(u)
end
