using MyFirstPackage
using Test

@testset "flip_direction_index" begin
    @test MyFirstPackage.flip_direction_index(D2Q9(), 1) == 9
end

