## Stingray: Next-Generation Spectral Timing
# Overview 

Stingray is a Julia library designed to perform time series analysis and related tasks on astronomical light curves. It supports a range of commonly used Fourier analysis techniques, as well as extensions for analyzing pulsar data, simulating datasets, and statistical modeling. Stingray is designed to be easy to extend and incorporate into data analysis workflows and pipelines.

# Features

## Current Capabilities

### 1. Data Handling and Simulation

- Good Time Interval (GTI) operations

### 2. Fourier Methods

- Basic Fourier transformations

### 3. Testing Framework

- Includes test cases for GTI and Fourier transformations
- Sample event list data for validation

## Future Plans

We welcome feature requests! If you need a particular tool that is not currently available or have a new method that could be useful in Stingray.jl, please contribute or open an issue.

### Upcoming Features

- **Power spectra calculations**
- **Extended power spectra and cross spectra support**
- **More comprehensive time series operations**
- **Advanced pulsar search methods**

# Installing Development Environment (For New Contributors)

For those of you wanting to contribute to the project, install the bleeding-edge development version from source. First, fork our repository on GitHub and clone the forked repository using:

```sh
$ git clone --recursive https://github.com/<your-github-username>/stingray.jl.git
```

Now, navigate to this folder and run the following command to add an upstream remote that’s linked to Stingray’s main repository. (This will be necessary when submitting PRs later.):

```sh
$ cd stingray.jl
$ git remote add upstream https://github.com/StingraySoftware/stingray.jl.git
```

# Installation
### Setting up in VS Code

1. Open Julia REPL
   - Press `Ctrl + Shift + P`
   - Search for `julia:start REPL` and select it
2. Navigate to the Cloned Repository in REPL:

```julia
cd("C:\\Users\\Stingray.jl")
```

3. Activate and Install Dependencies:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Test Suite

Please be sure to run the test suite before you use the package, and report anything you think might be bugs on our GitHub Issues page.

```julia
using Pkg
Pkg.test("Stingray")
```