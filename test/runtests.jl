using MyFirstPackage
using Test

@testset "MyFirstPackage.jl" begin
    include("lorenz.jl")
end
