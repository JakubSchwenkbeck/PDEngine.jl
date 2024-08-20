module PDEngine

# Import necessary packages
using Printf  # For formatted printing

# Include PDE solver files
include("heat_eq.jl")
include("wave_eq.jl")
include("navier_stokes_eq.jl")
include("poisson_eq.jl")

# Export functions for external use
export heat_eq_1d_fdm,heat_eq_1d_fem,
 wave_eq_1d_fdm,wave_eq_1d_fem
  navier_stokes_2d, 
  poisson_2d

# Wrapper functions or additional setup 
# initialize some parameters or configurations

# Example function - basic information about the package
function info()
    @info "PDEngine.jl: A Julia library for solving partial differential equations (PDEs)."
    @info "Available solvers:"
    @info "- heat_eq_1d: Solves the 1D heat equation."
    @info "- wave_eq_1d: Solves the 1D wave equation."
    @info "- navier_stokes_2d: Solves the 2D Navier-Stokes equations."
    @info "- poisson_2d: Solves the 2D Poisson's equation."
end

end # module PDEngine
