using Stingray
using Test
using FFTW, Distributions, Statistics, StatsBase, HDF5, FITSIO
using Logging, LinearAlgebra
using CFITSIO
using Random
using Plots

<<<<<<< HEAD
include("test_fourier.jl")
include("test_rebinning.jl")
=======
@testset "Fourier" begin
    include("test_fourier/test_fourier.jl")
    include("test_fourier/test_coherence.jl")
    include("test_fourier/test_norm.jl")
end
>>>>>>> upstream/unstable
@testset "GTI" begin
    include("test_gti.jl")
end
@testset "Eventlist" begin
    include("test_events.jl")
end
@testset verbose = true "Synthetic Events Tests" begin
    include("test_missionSupport.jl")
end
@testset "lightcurve" begin
    include("test_lightcurve.jl")
end
@testset "powerspectrum" begin
    include("test_powerspectrum.jl")
end
@testset "crossspectrum" begin
    include("test_crossspectrum.jl")
end

@testset "recipes" begin
    include("test_plotting/test_plots_recipes_lightcurve.jl")
    include("test_plotting/test_plots_recipes_gti.jl")
    include("test_plotting/test_plots_recipes_powerspectrum.jl")
    include("test_plotting/test_plots_recipes_crossspectrum.jl")
end
