```@meta
CurrentModule = Stingray
```

# Stingray

Documentation for [Stingray](https://github.com/matteobachetti/Stingray.jl).

```@index
```

```@autodocs
Modules = [Stingray]
```
# Stingray.jl

Stingray.jl is a spectral-timing software package for time series analysis of astronomical data, with a particular focus on X-ray astronomy.

## Installation

### Prerequisites

Before installing Stingray.jl, ensure you have Julia v1.7 or later installed:

1. Download Julia from the [official website](https://julialang.org/downloads/)
2. Follow the installation instructions for your operating system
3. Verify your installation by running `julia --version` in your terminal

### Installing from the Julia General Registry

To install the stable version of Stingray.jl:

```julia
using Pkg
Pkg.add("Stingray")
```

### Installing the Development Version

For the latest features and improvements, you can install directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url="https://github.com/matteobachetti/Stingray.jl")
```

Or switch to development mode if you plan to contribute:

```julia
using Pkg
Pkg.develop("Stingray")
```

### Verifying Installation

Verify your installation works correctly:

```julia
using Stingray
@info "Stingray.jl loaded successfully!"
```

## Testing

Stingray.jl includes a comprehensive test suite to ensure functionality. To run the tests:

```julia
using Pkg
Pkg.test("Stingray")
```

For contributors working on a local development version:

```julia
# Navigate to the package directory
cd(joinpath(Pkg.devdir(), "Stingray"))

# Run tests
using Pkg
Pkg.test()
```

### Writing New Tests

When adding features or fixing bugs, please include appropriate tests:

1. Tests should be placed in the `test/` directory
2. Use the `@testset` macro to organize related tests
3. Follow existing naming conventions (e.g., `test_[module_name].jl` or `test_[feature].jl`)
4. Run tests locally before submitting a pull request

Example test structure:

```julia
@testset "Module Name Tests" begin
    @test function_to_test(input) == expected_output
    @test_throws ErrorType function_to_test(invalid_input)
end
```

## Documentation

### Building Documentation Locally

To build and preview the documentation locally:

1. Install required documentation packages:

```julia
using Pkg
Pkg.add(["Documenter", "DocumenterTools"])
```

2. Navigate to the docs directory and build:

```julia
cd("docs")
julia --project=. make.jl
```

3. Open `docs/build/index.html` in your browser

### Writing Documentation

Documentation uses [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/). When adding new features:

1. Add docstrings to functions, types, and modules following Julia's documentation guidelines
2. Use `"""..."""` for multi-line docstrings
3. Include Examples section with working code examples
4. For complex features, add a dedicated page in `docs/src/`

Example docstring:

```julia
"""
    power_color(frequency, power; kwargs...)

Calculate two power colors from a power spectrum.

# Arguments
- `frequency::Vector{Number}`: The frequencies of the power spectrum
- `power::Vector{Number}`: The power at each frequency

# Keywords
- `power_err=nothing`: Optional error measurements for power
- `freq_edges=[1/256, 1/32, 0.25, 2.0, 16.0]`: Frequency edges for calculations
- `return_log=false`: Return the base-10 logarithm of results

# Returns
- `pc0::Number`: The first power color
- `pc0_err::Number`: Error on the first power color
- `pc1::Number`: The second power color
- `pc1_err::Number`: Error on the second power color

# Examples
```julia
freq = collect(0.01:0.01:20)
power = 1 ./ freq
pc0, pc0_err, pc1, pc1_err = power_color(freq, power)
```

See Heil et al. 2015, MNRAS, 448, 3348 for details.
"""
function power_color(frequency, power; kwargs...)
    # Implementation
end
```

## Example Notebooks

Example Jupyter notebooks demonstrating Stingray.jl functionality are available in the [notebooks repository](https://github.com/kashish2210/notebooks).

These examples cover:
- Basic time series analysis
- Power spectral calculations
- Cross-spectral techniques
- Specialized techniques for X-ray astronomy

## Contributing

We welcome contributions! Please check the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Implement your changes with tests and documentation
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to your branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Stingray.jl is a port of the Python [Stingray](https://github.com/StingraySoftware/stingray) package with optimizations for Julia.