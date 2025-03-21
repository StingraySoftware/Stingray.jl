using Test
using Stingray.Exceptions 

@testset "Exceptions Tests" begin
    err = StingrayError("This is error")
    
    @test err.msg == "This is error" 
    @test_throws StingrayError throw(StingrayError("Error test"))  
end
