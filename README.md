# Stingray.jl
Julia porting of Stingray

Stingray.jl is a Julia package for spectral timing analysis in X-ray astronomy. It provides functionality for power spectral density (PDS), cross-spectrum analysis, and variability estimation.

---

## 🚀 Installation

You can install **Stingray.jl** using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/StingraySoftware/Stingray.jl")
```

or, if you are working in a development environment:

```julia
Pkg.develop(url="https://github.com/StingraySoftware/Stingray.jl")
```

To ensure all dependencies are installed, run:

```julia
Pkg.instantiate()
```

---

## 📌 Usage

Once installed, you can start using **Stingray.jl** in your Julia environment:

```julia
using Stingray
```

---

## 🛠 Running Tests

To verify that the package is working correctly, navigate to the project directory and run:

```julia
using Pkg
Pkg.test("Stingray")
```

This will run all test cases located in the `test/` directory.

---

## 📖 Generating Documentation

To build the documentation, navigate to the `docs/` directory:

```sh
cd docs
```

Then run the following command in Julia:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("make.jl")
```

This will generate the documentation, which you can view in the `docs/build/` directory.

---

## 📂 Project Structure

```
Stingray.jl/
│── src/         # Main source code
│── test/        # Unit tests
│── docs/        # Documentation
│── Project.toml # Package dependencies
│── README.md    # This file
│── runtests.jl  # Test execution file
```

---

## Contributing

We welcome contributions! If you find an issue or want to add a feature, feel free to:

1. **Fork** the repository.
2. **Clone** your fork:  
   ```sh
   git clone https://github.com/your-username/Stingray.jl.git
   ```
3. **Create a feature branch**:
   ```sh
   git checkout -b new-feature
   ```
4. **Make changes and commit**:
   ```sh
   git commit -m "Added new feature"
   ```
5. **Push to your fork and create a PR**.

---

##  Important Instructions

- Ensure **Julia version >= 1.10**.
- Use **`Pkg.instantiate()`** before running tests or documentation.
- **FITS file support** requires `FITSIO.jl`.
- **Documentation uses a separate package environment** in `docs/`.
- Use **`Pkg.test("Stingray")`** instead of running individual test scripts.

---

### 🔗 Related Links

- 🏠 **GitHub Repository**: [Stingray.jl](https://github.com/StingraySoftware/Stingray)
