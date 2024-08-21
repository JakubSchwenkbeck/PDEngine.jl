module PDEngine


# Include PDE solver files
include("heat_eq.jl")
include("wave_eq.jl")
include("navier_stokes_eq.jl")
include("poisson_eq.jl")

# Export functions for external use
export heat_fdm, heat_fem,heat_nicolson,heat_spectral, wave_fdm, wave_fem,navier_stokes, poisson_fdm, poisson_fem 

# Wrapper functions or additional setup 
# initialize some parameters or configurations

#  basic information about the package
function info()
    @info "PDEngine.jl: A Julia library for solving partial differential equations (PDEs)."
    @info "Available solvers:"
    @info "- heat_fdm: Solves the 1D heat equation with the finite difference method"
    @info "- heat_fem: Solves the 1D heat equation with the finite element method"
    @info "- heat_nicolson: Solves the 1D heat equation with the crank nicolson method"
    @info "- wave_fdm: Solves the 1D wave equation with the finite difference method"
    @info "- wave_fem: Solves the 1D wave equation with the finite element method"
    @info "- navier_stokes: Solves the 2D Navier-Stokes equations."
    @info "- poisson_fdm: Solves the 2D Poisson's equation with the finite difference method"
    @info "- poisson_fem: Solves the 2D Poisson's equation with the finite element method"
end

end # module PDEngine
