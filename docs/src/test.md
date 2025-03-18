# Testing the Stingray Package

## Prerequisites
Before running the tests, ensure you have:
- **Julia 1.7 or later** installed.
- **Git** installed (for cloning the repository).
- All dependencies correctly installed.

## Installation and Setup
### Step 1: Clone the Stingray Repository
Open a terminal and run:
```sh
git clone https://github.com/your-username/Stingray.jl.git
cd Stingray.jl
```

### Step 2: Activate the Julia Environment
Inside the **Stingray** directory, launch the Julia REPL by running:
```sh
julia --project=.
```

### Step 3: Install Dependencies
Inside the Julia REPL, press `]` to enter **package mode**, then run:
```julia
instantiate
```
This will install all necessary dependencies.

## Running the Tests
### Option 1: Using Julia's Pkg Mode
1. **Open Julia** inside the project directory:
   ```sh
   julia --project=.
   ```
2. **Enter package mode** (press `]`) and run:
   ```julia
   test Stingray
   ```

### Option 2: Running `runtests.jl` Directly
Alternatively, you can run the test script manually:
```sh
julia --project=. test/runtests.jl
```

## Linux/macOS Commands Summary
For quick reference, here are all the necessary commands:
```sh
# Clone the repository
git clone https://github.com/your-username/Stingray.jl.git
cd Stingray.jl

# Open Julia with the project environment
julia --project=.

# Install dependencies (inside Julia REPL)
] instantiate

# Run tests (inside Julia REPL)
] test Stingray

# Or run tests directly
julia --project=. test/runtests.jl
```

## Troubleshooting
- If a package is missing, run:
  ```julia
  ] resolve
  ```
- If `runtests.jl` fails due to missing dependencies, ensure they are correctly installed:
  ```julia
  ] instantiate
  ```
- If tests fail, check for recent updates in the repository and run:
  ```sh
  git pull
  ```

## Final Notes
- This document assumes that `runtests.jl` includes tests for **all functions**.
- If you add new functions, ensure corresponding tests are included in `test/` before running `Pkg.test("Stingray")`.

