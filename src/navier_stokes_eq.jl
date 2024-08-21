using LinearAlgebra

"""
Solve the 2D Navier-Stokes equations using the finite difference method.
The Navier-Stokes equations describe the motion of viscous fluid substances.

Arguments:
- N: Number of grid points in each direction (excluding boundary points)
- ν: Kinematic viscosity
- T: Total simulation time
- Δx: Spatial step size
- Δt: Temporal step size

Returns:
- u: x-component of velocity field
- v: y-component of velocity field
- p: Pressure field
"""
function navier_stokes(N, ν, T, Δx, Δt)
    # Number of grid points (including boundaries)
    N_plus = N + 1
    
    # Initialize velocity and pressure fields
    u = zeros(N_plus, N_plus)   # x-component of velocity
    v = zeros(N_plus, N_plus)   # y-component of velocity
    p = zeros(N_plus, N_plus)   # Pressure field
    u_new = similar(u)
    v_new = similar(v)
    
    # Precompute constants
    r = ν * Δt / Δx^2
    
    # Helper functions for boundary conditions
    function apply_boundary_conditions!(u, v, p)
        # Apply boundary conditions for u, v, p
        u[1, :] .= 0      # u = 0 at y = 0 (no-slip condition)
        u[end, :] .= 0    # u = 0 at y = 1 (no-slip condition)
        v[:, 1] .= 0      # v = 0 at x = 0 (no-slip condition)
        v[:, end] .= 0    # v = 0 at x = 1 (no-slip condition)
        # Pressure boundary conditions are implicitly handled in the pressure correction step
    end
    
    function compute_pressure_correction!(u, v, p)
        # Compute the pressure correction using a Poisson equation
        b = zeros(N_plus, N_plus)
        
        # Build the right-hand side of the pressure Poisson equation
        for i in 2:N
            for j in 2:N
                b[i, j] = (u[i, j] - u[i, j-1]) / Δx + (v[i, j] - v[i-1, j]) / Δx
            end
        end
        
        # Solve the Poisson equation for pressure correction
        for _ in 1:1000  # Simple iterative solver (Jacobi method for demonstration)
            p_old = copy(p)
            for i in 2:N
                for j in 2:N
                    p[i, j] = (1/4) * (p_old[i+1, j] + p_old[i-1, j] + p_old[i, j+1] + p_old[i, j-1] - b[i, j] * Δx^2)
                end
            end
        end
    end
    
    function update_velocity!(u, v, p)
        # Update u and v fields using pressure correction
        for i in 2:N
            for j in 2:N
                u_new[i, j] = u[i, j] - (Δt / Δx) * (p[i, j] - p[i, j-1])
                v_new[i, j] = v[i, j] - (Δt / Δx) * (p[i, j] - p[i-1, j])
            end
        end
        # Swap references to avoid unnecessary copying
        u, u_new = u_new, u
        v, v_new = v_new, v
    end

    # Time-stepping loop
    for _ in Δt:Δt:T
        # Compute intermediate velocities
        u_star = copy(u)
        v_star = copy(v)
        
        # Update velocities using the advection and diffusion terms
        for i in 2:N
            for j in 2:N
                u_star[i, j] = u[i, j] + Δt * (
                    - (u[i, j] * (u[i, j] - u[i, j-1]) / Δx + v[i, j] * (u[i, j] - u[i-1, j]) / Δx) +
                    ν * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / Δx^2 + (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / Δx^2)
                )
                
                v_star[i, j] = v[i, j] + Δt * (
                    - (u[i, j] * (v[i, j] - v[i, j-1]) / Δx + v[i, j] * (v[i, j] - v[i-1, j]) / Δx) +
                    ν * ((v[i+1, j] - 2 * v[i, j] + v[i-1, j]) / Δx^2 + (v[i, j+1] - 2 * v[i, j] + v[i, j-1]) / Δx^2)
                )
            end
        end
        
        # Apply boundary conditions
        apply_boundary_conditions!(u_star, v_star, p)
        
        # Compute pressure correction
        compute_pressure_correction!(u_star, v_star, p)
        
        # Update velocities using the pressure correction
        update_velocity!(u_star, v_star, p)
    end
    
    return u, v, p
end

# Example usage
N = 50
ν = 0.01
T = 1.0
Δx = 1.0 / N
Δt = 0.001

u, v, p = navier_stokes_2d(N, ν, T, Δx, Δt)
