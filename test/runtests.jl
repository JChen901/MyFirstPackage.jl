using MyFirstPackage
using Test

@testset "lorenz.jl" begin
    include("lorenz.jl")
end

@testset "fluid.jl" begin
    include("fluid.jl")
end