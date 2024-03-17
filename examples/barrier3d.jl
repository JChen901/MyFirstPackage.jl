using MyFirstPackage
using LinearAlgebra
using Makie: RGBA

lb = example_d3q19()

initial_state = collect(lb.grid[1].density)



points = Observable(Point3f[]) # Signal that can be used to update plots efficiently
colors = Observable(Int[])

set_theme!(theme_black())

fig, ax, l = lines(points, color = colors,
	colormap = :inferno, transparency = true,
	axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
		viewmode = :fit, limits = (0, 40, 0, 50, 0, 50)))

record(fig, "lorenz.mp4", 1:200) do frame
	if frame % 10 == 0
		step!(lb)
		@show frame
		for ci in CartesianIndices(lb.grid)
			p = lb.grid[ci]
			if lb.barrier[ci]
				push!(points[], ci.I)
				push!(colors[], 500)
			end
			diff = norm(collect(p.density) - initial_state)
			if diff > 0.05
				push!(points[], ci.I)
				push!(colors[], round(diff * 1000))
			end
		end
	end
	ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120) # set the view angle of the axis
	notify(points)
	notify(colors) # tell points and colors that their value has been updated
	l.colorrange = (0, frame) # update plot attribute directly
end
