using Test
using MyFirstPackage

@testset "directions" begin
    for i in 1:19
        @test directions(D3Q19())[i]+directions(D3Q19())[MyFirstPackage.flip_direction_index(D3Q19(), i)] == Point(0, 0, 0)
    end
end

@testset "LatticeBoltzmann amd stream" begin
    h = 80
    l = 100
    w = 200
	u0 = Point(0.0, 0.1,0.2)
	rho = equilibrium_density(D3Q19(), 1.0, u0)
	rgrid = fill(rho, (h,l, w))
	# Initialize barriers:
	barrier = falses(h, l,w)  # True wherever there's a barrier
	lb = LatticeBoltzmann3d(D3Q19(), rgrid, barrier)
    step!(lb)
end