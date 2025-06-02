"""
Dictionary of simple conversion functions for different missions.

This dictionary provides PI (Pulse Invariant) to energy conversion functions
for various X-ray astronomy missions. Each function takes a PI channel value
and returns the corresponding energy in keV.

Supported missions:
- NuSTAR: Nuclear Spectroscopic Telescope Array
- XMM: X-ray Multi-Mirror Mission  
- NICER: Neutron star Interior Composition Explorer
- IXPE: Imaging X-ray Polarimetry Explorer
- AXAF/Chandra: Advanced X-ray Astrophysics Facility
- XTE/RXTE: Rossi X-ray Timing Explorer
"""
const SIMPLE_CALIBRATION_FUNCS = Dict{String, Function}(
    "nustar" => (pi) -> pi * 0.04 + 1.62,
    "xmm" => (pi) -> pi * 0.001,
    "nicer" => (pi) -> pi * 0.01,
    "ixpe" => (pi) -> pi / 375 * 15,
    "axaf" => (pi) -> (pi - 1) * 14.6e-3,  # Chandra/AXAF
    "chandra" => (pi) -> (pi - 1) * 14.6e-3,  # Explicit chandra entry
    "xte" => (pi) -> pi * 0.025  # RXTE/XTE
)

"""
Abstract type for mission-specific calibration and interpretation.

This serves as the base type for all mission support implementations,
allowing for extensibility and type safety in mission-specific operations.
"""
abstract type AbstractMissionSupport end

"""
    MissionSupport{T} <: AbstractMissionSupport

Structure containing mission-specific calibration and interpretation information.

This structure encapsulates all the necessary information for handling
data from a specific X-ray astronomy mission, including calibration
functions, energy column alternatives, and GTI extension preferences.

# Fields
- `name::String`: Mission name (normalized to lowercase)
- `instrument::Union{String, Nothing}`: Instrument identifier
- `epoch::Union{T, Nothing}`: Observation epoch in MJD (for time-dependent calibrations)
- `calibration_func::Function`: PI to energy conversion function
- `interpretation_func::Union{Function, Nothing}`: Mission-specific FITS interpretation function
- `energy_alternatives::Vector{String}`: Preferred energy column names in order of preference
- `gti_extensions::Vector{String}`: GTI extension names in order of preference

# Type Parameters
- `T`: Type of the epoch parameter (typically Float64)
"""
struct MissionSupport{T} <: AbstractMissionSupport
    name::String
    instrument::Union{String, Nothing}
    epoch::Union{T, Nothing}
    calibration_func::Function
    interpretation_func::Union{Function, Nothing}
    energy_alternatives::Vector{String}
    gti_extensions::Vector{String}
end

"""
    get_mission_support(mission::String, instrument=nothing, epoch=nothing) -> MissionSupport

Create mission support object with mission-specific parameters.

This function creates a MissionSupport object containing all the necessary
information for processing data from a specified X-ray astronomy mission.
It handles mission aliases (e.g., Chandra/AXAF) and provides appropriate
defaults for each mission.

# Arguments
- `mission::String`: Mission name (case-insensitive)
- `instrument::Union{String, Nothing}=nothing`: Instrument identifier
- `epoch::Union{Float64, Nothing}=nothing`: Observation epoch in MJD

# Returns
- `MissionSupport{Float64}`: Mission support object

# Throws
- `ArgumentError`: If mission name is empty

# Examples
```julia
# Basic usage
ms = get_mission_support("nustar")

# With instrument specification
ms = get_mission_support("nustar", "FPM_A")

# With epoch for time-dependent calibrations
ms = get_mission_support("xte", "PCA", 50000.0)
```
"""
function get_mission_support(mission::String, 
                           instrument::Union{String, Nothing}=nothing,
                           epoch::Union{Float64, Nothing}=nothing)
    
    # Check for empty mission string
    if isempty(mission)
        throw(ArgumentError("Mission name cannot be empty"))
    end
    
    mission_lower = lowercase(mission)
    
    # Handle chandra/axaf aliases - normalize to chandra
    if mission_lower in ["chandra", "axaf"]
        mission_lower = "chandra"
    end
    
    calib_func = if haskey(SIMPLE_CALIBRATION_FUNCS, mission_lower)
        SIMPLE_CALIBRATION_FUNCS[mission_lower]
    else
        @warn "Mission $mission not recognized, using identity function"
        identity
    end
    
    # Mission-specific energy alternatives (order matters!)
    energy_alts = if mission_lower in ["chandra", "axaf"]
        ["ENERGY", "PI", "PHA"]  # Chandra usually has ENERGY column
    elseif mission_lower == "xte"
        ["PHA", "PI", "ENERGY"]
    elseif mission_lower == "nustar"
        ["PI", "ENERGY", "PHA"]
    else
        ["ENERGY", "PI", "PHA"]
    end
    
    # Mission-specific GTI extensions
    gti_exts = if mission_lower == "xmm"
        ["GTI", "GTI0", "STDGTI"]
    elseif mission_lower in ["chandra", "axaf"]
        ["GTI", "GTI0", "GTI1", "GTI2", "GTI3"]
    else
        ["GTI", "STDGTI"]
    end
    
    MissionSupport{Float64}(mission_lower, instrument, epoch, calib_func, nothing, energy_alts, gti_exts)
end

"""
    apply_calibration(mission_support::MissionSupport, pi_channels::AbstractArray) -> Vector{Float64}

Apply calibration function to PI channels.

Converts PI (Pulse Invariant) channel values to energies in keV using
the mission-specific calibration function stored in the MissionSupport object.

# Arguments
- `mission_support::MissionSupport`: Mission support object containing calibration function
- `pi_channels::AbstractArray{T}`: Array of PI channel values

# Returns
- `Vector{Float64}`: Array of energy values in keV

# Examples
```julia
ms = get_mission_support("nustar")
pi_values = [100, 500, 1000]
energies = apply_calibration(ms, pi_values)
```
"""
function apply_calibration(mission_support::MissionSupport, pi_channels::AbstractArray{T}) where T
    if isempty(pi_channels)
        return similar(pi_channels, Float64)
    end
    return mission_support.calibration_func.(pi_channels)
end

"""
    patch_mission_info(info::Dict{String,Any}, mission=nothing) -> Dict{String,Any}

Apply mission-specific patches to header information.

This function applies mission-specific modifications to FITS header information
to handle mission-specific quirks and conventions. It's based on the Python
implementation in Stingray's mission interpretation module.

# Arguments
- `info::Dict{String,Any}`: Dictionary containing header information
- `mission::Union{String,Nothing}=nothing`: Mission name

# Returns
- `Dict{String,Any}`: Patched header information dictionary

# Examples
```julia
info = Dict("gti" => "STDGTI", "ecol" => "PHA")
patched = patch_mission_info(info, "xmm")  # Adds GTI0 to gti field
```
"""
function patch_mission_info(info::Dict{String,Any}, mission::Union{String,Nothing}=nothing)
    if isnothing(mission)
        return info
    end
    
    mission_lower = lowercase(mission)
    patched_info = copy(info)
    
    # Normalize chandra/axaf
    if mission_lower in ["chandra", "axaf"]
        mission_lower = "chandra"
    end
    
    if mission_lower == "xmm" && haskey(patched_info, "gti")
        patched_info["gti"] = string(patched_info["gti"], ",GTI0")
    elseif mission_lower == "xte" && haskey(patched_info, "ecol")
        patched_info["ecol"] = "PHA"
        patched_info["ccol"] = "PCUID"
    elseif mission_lower == "chandra"
        # Chandra-specific patches
        if haskey(patched_info, "DETNAM")
            patched_info["detector"] = patched_info["DETNAM"]
        end
        # Add Chandra-specific time reference if needed
        if haskey(patched_info, "TIMESYS")
            patched_info["time_system"] = patched_info["TIMESYS"]
        end
    end
    
    return patched_info
end
function interpret_fits_data!(f::FITS, mission_support::MissionSupport)
    # Placeholder for mission-specific interpretation
    # This would contain mission-specific FITS handling logic
    return nothing
end