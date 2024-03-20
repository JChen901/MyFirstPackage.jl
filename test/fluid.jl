using Test
using MyFirstPackage

@testset "directions" begin
    for i in 1:9
        @test directions(D2Q9())[i]+directions(D2Q9())[MyFirstPackage.flip_direction_index(D2Q9(), i)] == Point(0, 0)
    end
end

@testset "momentum" begin
    @test momentum(D2Q9(), MyFirstPackage.Cell((fill(1,9)...,))) == Point(0.0, 0.0)
    @test momentum(D2Q9(), MyFirstPackage.Cell((1,fill(0,8)...,))) == Point(1.0, 1.0)
end

@testset "LatticeBoltzmann amd stream" begin
    height = 80
    width = 200
	u0 = Point(0.0, 0.1)
	rho = equilibrium_density(D2Q9(), 1.0, u0)
	rgrid = fill(rho, height, width)
	# Initialize barriers:
	barrier = falses(height, width)  # True wherever there's a barrier
	mid = div(height, 2)
	barrier[mid-8:mid+8, div(height, 2)] .= true              # simple linear barrier
	lb = LatticeBoltzmann(D2Q9(), rgrid, barrier)
    MyFirstPackage.stream!(lb.config, lb.grid, lb.gridcache, lb.barrier)
    lb.grid .= MyFirstPackage.collide.(Ref(lb.config), lb.grid)
end
