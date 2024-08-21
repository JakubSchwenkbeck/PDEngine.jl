module PDEngine


# Include PDE solver files
include("heat_eq.jl")
include("wave_eq.jl")
include("navier_stokes_eq.jl")
include("poisson_eq.jl")

# Export functions for external use
export heat_eq_1d_fdm, heat_eq_1d_fem,heat_eq_1d_crank_nicolson, wave_eq_1d_fdm, wave_eq_1d_fem,navier_stokes_2d, poisson_2d_fdm, poisson_2d_fem 

# Wrapper functions or additional setup 
# initialize some parameters or configurations

#  basic information about the package
function info()
    @info "PDEngine.jl: A Julia library for solving partial differential equations (PDEs)."
    @info "Available solvers:"
    @info "- heat_eq_1d_fdm: Solves the 1D heat equation with the finite difference method"
    @info "- heat_eq_1d_fem: Solves the 1D heat equation with the finite element method"
    @info "- wave_eq_1d_fdm: Solves the 1D wave equation with the finite difference method"
    @info "- wave_eq_1d_fem: Solves the 1D wave equation with the finite element method"
    @info "- navier_stokes_2d: Solves the 2D Navier-Stokes equations."
    @info "- poisson_2d_fdm: Solves the 2D Poisson's equation with the finite difference method"
    @info "- poisson_2d_fem: Solves the 2D Poisson's equation with the finite element method"
end

end # module PDEngine
