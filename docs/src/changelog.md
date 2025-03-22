---
# Update on Dependency to Use Julia Version > 1.10  

## Changelog  

### Summary  
This update refactors the **Stingray.jl** package to ensure compatibility with **Julia versions 1.10 and above** while improving functionality, test coverage, and dependency management.  

### Changes Implemented  

1. **Dependency Updates**  
   - Updated `Project.toml` to enforce **Julia 1.10+** compatibility.  
   - Removed unnecessary dependencies and streamlined package imports.  
   - Ensured **`ResumableFunctions`** is correctly included in both `[deps]` and `[extras]`.  

2. **Improved Functionality**  
   - Introduced **`sum_if_not_none_or_initialize`** to prevent array interface issues:
     ```julia
     # This function ensures that the sum operation maintains the array structure
     # instead of returning a scalar sum, avoiding potential type mismatches.
     ```
   - Enhanced performance of **cross-spectrum** and **power spectral density** calculations.  

3. **Testing Enhancements**  
   - Consolidated all test cases into **`runtests.jl`** under the `test/` directory.  
   - Added missing **normalization tests** (`test_norm`).  
   - Included test commands for **Linux and macOS** in `test.md`.  

### References  
- Dependency Fixes: [Stingray.jl Issue #17](https://github.com/StingraySoftware/Stingray.jl/issues/17)  
- Followed by: [Stingray.jl Issue #16](https://github.com/StingraySoftware/Stingray.jl/issues/16)  

---

