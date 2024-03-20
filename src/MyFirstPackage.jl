module MyFirstPackage
# import packages
using LinearAlgebra

# export interfaces
export Lorenz, integrate_step
export Point, Point2D, Point3D
export RungeKutta, Euclidean
export D2Q9, LatticeBoltzmann, step!, equilibrium_density, momentum, curl, example_d2q9, density, directions 

export D3Q19,LatticeBoltzmann3d,example_d3q19

include("lorenz.jl")
include("fluid.jl")
include("fluid3d.jl")
end
