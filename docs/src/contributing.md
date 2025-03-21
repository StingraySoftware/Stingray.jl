# Installation Guide for Stingray.jl
Documentation for [Stingray.jl](https://github.com/StingraySoftware/Stingray.jl).

## installing Julia

1. Download Julia from the official site: [https://julialang.org/downloads/](https://julialang.org/downloads/)
2. Install it following your OS instructions.
3. Add Julia to your system's PATH to summon it from any terminal like magic.

---

## Installing Stingray.jl

Launch the Julia REPL, press `]` to enter the **package manager mode**, and run:

```julia
using Pkg
Pkg.add("Stingray")
```

This will fetch Stingray from the heavens (and install dependencies).

---

## Running Tests Locally

To test the waters and ensure the installation vibes are right, run this in the **main directory**:

```julia
using Pkg
Pkg.test("Stingray")
```

If you get a dependency error with Stingray:

```julia
Pkg.add(url="https://github.com/StingraySoftware/Stingray.jl")
Pkg.develop("Stingray")
```

---

## Cloning the Development Version

```julia
Pkg.add(url="https://github.com/matteobachetti/Stingray.jl")
Pkg.develop("Stingray")
```

Or if you already cloned it manually:

```julia
using Pkg
Pkg.develop(path="path/to/your/local/Stingray.jl")
```

---

## Building Documentation (the doc-gen ritual)

1. Navigate to the `docs/` folder:
   ```bash
   cd docs
   ```

2. Activate the project:
   ```julia
   using Pkg
   Pkg.activate(".")
   ```

3. Build the docs:
   ```julia
   include("make.jl")
   ```


##  Checking Your Installation

```julia
using Stingray
@info "Stingray.jl loaded successfully!"
```

---

## error handling

- If you’re juggling environments, run this before installing anything:
  ```julia
  Pkg.activate("path/to/your/project")
  ```

- If you want to **remove and reinstall** for a clean slate:
  ```julia
  Pkg.rm("Stingray")
  Pkg.add("Stingray")
  ```

- To update the package:
  ```julia
  Pkg.update("Stingray")
  ```

- To see all dependencies and versions:
  ```julia
  Pkg.status()
  ```

---

**Julia porting of Stingray - Under heavy development, be ready to help debugging it ;)**
