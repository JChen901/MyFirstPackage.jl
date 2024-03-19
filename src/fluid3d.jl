"""
	D3Q19 <: AbstractLBConfig{3, 19}

A lattice Boltzmann configuration for 3D, 19-velocity model.
"""
struct D3Q19 <: AbstractLBConfig{3, 19} end

directions(::D3Q19) = (
	Point(1, 0, -1),
	Point(0, 1, 1),
	Point(0, 1, -1),
	Point(1, 0, 1),
	Point(1, -1, 0),
	Point(1, 1, 0),
	Point(0, 0, 1),
	Point(0, 1, 0),
	Point(1, 0, 0),
	Point(0, 0, 0),
	Point(-1, 0, 0),
	Point(0, -1, 0),
	Point(0, 0, -1),
	Point(-1, -1, 0),
	Point(-1, 1, 0),
	Point(-1, 0, -1),
	Point(0, -1, 1),
	Point(0, -1, -1),
	Point(-1, 0, 1),
)

# directions[k] is the opposite of directions[flip_direction_index(k)
function flip_direction_index(::D3Q19, i::Int)
	return 20 - i
end

weights(::D3Q19) = (1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 18, 1 / 18, 1 / 18, 1 / 3, 1 / 18, 1 / 18, 1 / 18, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36)

"""
	LatticeBoltzmann{D, N, T, CFG, MT, BT}

A lattice Boltzmann simulation with D dimensions, N velocities, and lattice configuration CFG.
"""
struct LatticeBoltzmann3d{D, N, T, CFG <: AbstractLBConfig{D, N}, MT <: AbstractArray{Cell{N, T}}, BT <: AbstractArray{Bool}}
	config::CFG # lattice configuration
	grid::MT    # density of the fluid
	gridcache::MT # cache for the density of the fluid
	barrier::BT # barrier configuration
end

function LatticeBoltzmann3d(config::AbstractLBConfig{D, N}, grid::AbstractArray{<:Cell}, barrier::AbstractArray{Bool}) where {D, N}
	@assert size(grid) == size(barrier)
	return LatticeBoltzmann3d(config, grid, similar(grid), barrier)
end


function stream!(
	lb::AbstractLBConfig{3, N},  # lattice configuration for 3D
	newgrid::AbstractArray{D}, # the updated grid (3D)
	grid::AbstractArray{D}, # the original grid (3D)
	barrier::AbstractArray{Bool}, # the barrier configuration (3D)
) where {N, T, D <: Cell{N, T}}
	ds = directions(lb)  # Ensure this provides the D3Q19 velocity vectors
	@inbounds for ci in CartesianIndices(newgrid)
		i, j, k = ci.I  # Now we have three indices for 3D
		newgrid[ci] = Cell(ntuple(N) do d  # Collect the densities
			ei = ds[d]
			l, m, n = size(grid)  # Dimensions for 3D grid
			i2, j2, k2 = mod1(i - ei[1], l), mod1(j - ei[2], m), mod1(k - ei[3], n)
			if barrier[i2, j2, k2]
				# if the cell is a barrier, the fluid flows back
				density(grid[i, j, k], flip_direction_index(lb, d))
			else
				# otherwise, the fluid flows to the neighboring cell
				density(grid[i2, j2, k2], d)
			end
		end)
	end
end

function step!(lb::LatticeBoltzmann3d)
	copyto!(lb.gridcache, lb.grid)
	stream!(lb.config, lb.grid, lb.gridcache, lb.barrier)
	lb.grid .= collide.(Ref(lb.config), lb.grid)
	return lb
end

function example_d3q19(;
	height = 40, width = 50, length = 50,
	u0 = Point(0.1, 0.0, 0.0)) # initial and in-flow speed
	# Initialize all the arrays to steady rightward flow:
	rho = equilibrium_density(D3Q19(), 1.0, u0)
	rgrid = fill(rho, (height, width, length))

	# Initialize barriers:
	barrier = falses(height, width, length)  # True wherever there's a barrier
	mid = div(width, 2)
	barrier[div(height, 2), mid-8:mid+8, mid-8:mid+8] .= true              # simple linear barrier

	return LatticeBoltzmann3d(D3Q19(), rgrid, barrier)
end

function curl(u::Matrix{Point3D{T}}) where T
	return map(CartesianIndices(u)) do ci
		i, j, k = ci.I
		m, n, l = size(u)
		uypx = u[mod1(i + 1, m), j,k][2] - u[mod1(i - 1, m), j,k][2]
		uxpy = u[i, mod1(j + 1, n),k][1] - u[i, mod1(j - 1, n),k][1]

		uypz = u[i, j, mod1(k + 1, l)][2] - u[i, j, mod1(k - 1, l)][2]
		uzpy = u[i, mod1(j + 1, n),k][3] - u[i, mod1(j - 1, n),k][3]

		uzpx = u[mod1(i + 1, m), j,k][3] - u[mod1(i - 1, m), j,k][3]
		uxpz = u[i, j, mod1(k + 1, l)][1] - u[i, j, mod1(k - 1, l)][1]
		return uzpy - uypz, uxpz - uzpx, uypx - uxpy
	end
end
