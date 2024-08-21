"""
Solve 2D Poisson's equation using the finite difference method.
Poisson's equation is given by: ∇²u = -f
where f is a source term.
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
Poisson's equation is given by: ∇²u = -f
where f is a source term.

Arguments:
- N: Number of grid points along each dimension (excluding boundaries)
- f: A function representing the source term
- Δx: Spatial step size

Returns:
- u: The final solution field
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