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
