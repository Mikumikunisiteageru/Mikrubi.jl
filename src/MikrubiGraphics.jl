module MikrubiGraphics

using PyPlot
using GeoArrays
using Shapefile
using StatsBase

export showlayer, showfield, showctpixels, showshptable

"""
	showlayer(layer; f=identity, kwargs...)

Shows a layer. Keyword argument `f = identity` is a function acted separately
on every element. A possible alternative is `f = x -> x ^ 0.4`. 
"""
function showlayer(layer; f=identity, kwargs...)
	cc = coords(layer)
	x = first.(cc)
	y = last.(cc)
	A = (x -> ismissing(x) ? NaN : x).(layer.A[:, :, 1])
	nna = findall(.!isnan.(A))
	lim(z) = [1.05 -0.05; -0.05 1.05] * 
		[minimum(vcat(z[nna], z[nna .+ [CartesianIndex(1,1)]])), 
		 maximum(vcat(z[nna], z[nna .+ [CartesianIndex(1,1)]]))]
	pcolor(x, y, f.(A); kwargs...)
	xlim(lim(x))
	ylim(lim(y))
end

"""
	showlayer(layer; f=identity, kwargs...)
	showfield(field, layer; f=tiedrank, kwargs...)

Shows geographic information and environmental information of a Mikrubi model.
The three principal components are reflexed in red, green, and blue. Keyword
argument `f = tiedrank` is a function acted on columns of `field.vars` as a 
whole. A possible alternative is `f = identity`.
"""
function showfield(field; f=tiedrank, kwargs...)
	u(a) = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
	v(a) = u(f(a))
	r = field.dvar >= 1 ? v(field.vars[:, 1]) : fill(0.5, field.npixel)
	g = field.dvar >= 2 ? v(field.vars[:, 2]) : fill(0.5, field.npixel)
	b = field.dvar >= 3 ? v(field.vars[:, 3]) : fill(0.5, field.npixel)
	@assert size(field.locs, 2) >= 2
	x = field.locs[:, 1]
	y = field.locs[:, 2]
	scatter(x, y, c=collect(zip(r, g, b)), s=2; kwargs...)
end
function showfield(field, layer; f=tiedrank, kwargs...)
	u(a) = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
	v(a) = u(f(a))
	r = field.dvar >= 1 ? v(field.vars[:, 1]) : fill(0.5, field.npixel)
	g = field.dvar >= 2 ? v(field.vars[:, 2]) : fill(0.5, field.npixel)
	b = field.dvar >= 3 ? v(field.vars[:, 3]) : fill(0.5, field.npixel)
	@assert layer.f.linear[1,2] == layer.f.linear[2,1] == 0
	xb, yb, _ = size(layer)
	matrix = fill(1., yb, xb, 3)
	for i = 1:field.npixel
		# x, y = indices(layer, field.locs[i, :])
		x, y = ceil.(Int, inv(layer.f)(field.locs[i, :]))
		matrix[y, x, :] .= r[i], g[i], b[i]
	end
	x1, y1 = coords(layer, [1, 1])
	x2, y2 = coords(layer, [xb+1, yb+1])
	imshow(matrix, extent=(x1,x2, y2,y1), zorder=-1, aspect="auto"; kwargs...)
	val = minimum(matrix, dims=3) .< 1
	xval = any(val, dims=1)[:]
	yval = any(val, dims=2)[:]
	x3, y3 = coords(layer, [findfirst(xval)-1, findfirst(yval)-1])
	x4, y4 = coords(layer, [findlast(xval), findlast(yval)])
	lim(z3, z4) = [1.05 -0.05; -0.05 1.05] * [min(z3, z4), max(z3, z4)]
	xlim(lim(x3, x4))
	ylim(lim(y3, y4))
end

pseudorand(x, y) = x % y / y
function hashcolor(x; salt=20, gmin=0.5)
	r0, g0, b0 = pseudorand.(hash((salt, x)), [39, 139, 239])
	r, g, b = @. 1 - (1 - (r0, g0, b0)) * (1 - gmin)
end

"""
	showctpixels(ctpixels, layer; salt=20, kwargs...)

Shows a `Mikrubi.CtPixels`. Every county is assigned a hash color (influenced
by a fixed `salt` value also), and every pixel has the composite color from all
counties assigned to it. Empty cells are depicted white.
"""
function showctpixels(ctpixels, layer; salt=20, kwargs...)
	@assert layer.f.linear[1,2] == layer.f.linear[2,1] == 0
	xb, yb, _ = size(layer)
	matrix = fill(1., yb, xb, 3)
	for ((x, y), id) = ctpixels
		matrix[y, x, :] .*= hashcolor(id, salt=salt)
	end
	x1, y1 = coords(layer, [1, 1])
	x2, y2 = coords(layer, [xb+1, yb+1])
	imshow(matrix, extent=(x1,x2, y2,y1), zorder=-1, aspect="auto"; kwargs...)
	val = minimum(matrix, dims=3) .< 1
	xval = any(val, dims=1)[:]
	yval = any(val, dims=2)[:]
	x3, y3 = coords(layer, [findfirst(xval)-1, findfirst(yval)-1])
	x4, y4 = coords(layer, [findlast(xval), findlast(yval)])
	lim(z3, z4) = [1.05 -0.05; -0.05 1.05] * [min(z3, z4), max(z3, z4)]
	xlim(lim(x3, x4))
	ylim(lim(y3, y4))
end

function polygonline(geom)
	line = Tuple{Float64, Float64}[]
	parts = geom.parts .+ 1
	for i = 1:length(geom.points)
		pt = geom.points[i]
		i in parts && push!(line, (NaN, NaN))
		push!(line, (pt.x, pt.y))
	end
	line
end
function shapelines(shptable)
	geoms = Shapefile.shapes(shptable)
	lines = Tuple{Float64, Float64}[]
	for geom = geoms
		append!(lines, polygonline(geom))
	end
	n = length(lines)
	flag = trues(n)
	set = Set([((NaN, NaN), (NaN, NaN))])
	for i = 2:n
		if (lines[i-1], lines[i]) in set
			flag[i] = false
			continue
		end
		push!(set, (lines[i-1], lines[i]))
		push!(set, (lines[i], lines[i-1]))
	end
	line = deepcopy(lines)
	line[.!flag] .= [(NaN, NaN)]
	shortline = Tuple{Float64, Float64}[]
	lastnan = true
	for i = 2:n
		isnan(line[i][1]) && lastnan && continue
		push!(shortline, line[i])
		lastnan = isnan(line[i][1])
	end
	shortline
end

"""
	showshptable(shptable; kwargs...)

Shows lines from polygons in `shptable`. Identical segments are reduced as one.
"""
function showshptable(shptable; kwargs...)
	shortline = shapelines(shptable)
	plot(first.(shortline), last.(shortline), "-k", lw=0.8; kwargs...)
end

end # module
