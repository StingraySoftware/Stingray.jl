using Test
using Pkg
using Plots

gr()

function get_package_version(pkg_name)
    for (_, pkg) in Pkg.dependencies()
        if pkg.name == pkg_name
            return pkg.version
        end
    end
    return "Not Found"
end

function show_tested_versions()
    tested_versions = Dict(
        "Julia" => VERSION,
        "Plots" => get_package_version("Plots")
    )

    println("Tested Versions:")
    for (pkg, ver) in tested_versions
        println("  $pkg: $ver")
    end
end

function runtests()
    show_tested_versions()
    @testset "Stingray.jl Tests" begin
        include(joinpath(@__DIR__, "..", "test", "test_fourier.jl"))
        include(joinpath(@__DIR__, "..", "test", "test_gti.jl"))
    end
end


runtests()
