using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging
using Pkg
using InteractiveUtils
# this function is for temporary support we can remove after the update of FFTW.jl
# to use the system FFTW library on macOS
function configure_fftw_for_macos()
    @info "Configuring FFTW for macOS..."
    fftw_env_vars = ["FFTW_LIBRARY_PATH", "FFTW_LIBRARY", "FFTW_FORCE_SYSTEM_LIB"]
    for var in fftw_env_vars
        if haskey(ENV, var)
            @info "Removing existing FFTW environment variable: $var"
            delete!(ENV, var)
        end
    end
    @info "Setting up FFTW with FFTW_jll provider..."
    FFTW.set_provider!("fftw")
    Pkg.build("FFTW")
end

if Sys.isapple()
    @info "macOS detected, applying FFTW configurations..."
    configure_fftw_for_macos()
end
include("test_fourier.jl")
include("test_gti.jl")
include("test_events.jl")
