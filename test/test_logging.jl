using Test
using Stingray

@testset "Logging Tests" begin
    @info "Testing INFO level log"
    @warn "Testing WARNING level log"
    @error "Testing ERROR level log"
    @test true  # Just a dummy test to ensure script runs
end
