# src/pyplot.jl

export showlayer, showfield, showctpixels, showshptable

function geom2mat(geom)
	GI.isgeometry(geom) || error("`geom` is not a geometry!")
	trait = GI.geomtrait(geom)
	if isa(trait, GI.PolygonTrait)
		return geom2mat_polygon(geom)
	elseif isa(trait, GI.MultiPolygonTrait)
		return geom2mat_multipolygon(geom)
	else
		error("Trait of `geom` not supported!")
	end
end

function geom2mat_polygon(polygon)
	parts = Int[]
	points = NTuple{2, Float64}[]
	for i = 1 : GI.nring(polygon)
		ring = GI.getring(polygon, i)
		push!(parts, length(points))
		for j = 1 : GI.npoint(ring)
			point = GI.getpoint(ring, j)
			push!(points, (GI.x(point), GI.y(point)))
		end
	end
	n = length(points)
	mat = Matrix{Float64}(undef, 2, n)
	for i = 1:n
		mat[:, i] .= points[i]
	end
	return parts, mat
end

function geom2mat_multipolygon(multipolygon)
	parts = Int[]
	points = NTuple{2, Float64}[]
	for k = 1 : GI.npolygon(multipolygon)
		polygon = GI.getpolygon(multipolygon, k)
		for i = 1 : GI.nring(polygon)
			ring = GI.getring(polygon, i)
			push!(parts, length(points))
			for j = 1 : GI.npoint(ring)
				point = GI.getpoint(ring, j)
				push!(points, (GI.x(point), GI.y(point)))
			end
		end
	end
	n = length(points)
	mat = Matrix{Float64}(undef, 2, n)
	for i = 1:n
		mat[:, i] .= points[i][1:2]
	end
	return parts, mat
end

function coords(grid::Raster)
	x, y = dims(grid)[1:2]
	xseq = dimseq(x)
	yseq = dimseq(y)
	# vcat.(xseq, yseq')
	xx = repeat(xseq, 1, length(yseq))
	yy = repeat(yseq', length(xseq), 1)
	return xx, yy
end

"""
	showlayer(layer; ax=PyPlot.gca(), f=identity, kwargs...)

Show a layer. Keyword argument `f = identity` is a function acted separately
on every element. A possible alternative is `f = x -> x ^ 0.4`. 
"""
function showlayer(layer; ax=PyPlot.gca(), f=identity, kwargs...)
	xx, yy = coords(layer)
	Anan = replace_missing(layer, missingval=eltype(layer)(NaN))
	A = Array(Anan)[:, :, 1]
	nna = findall(.!isnan.(A))
	lim(z) = [1.05 -0.05; -0.05 1.05] * 
		[minimum(vcat(z[nna], z[nna .+ [CartesianIndex(1,1)]])), 
		 maximum(vcat(z[nna], z[nna .+ [CartesianIndex(1,1)]]))]
	handle = ax.pcolor(xx, yy, f.(A); kwargs...)
	ax.set_xlim(lim(xx))
	ax.set_ylim(lim(yy))
	return handle
end

"""
	showfield(layer; ax=PyPlot.gca(), f=identity, kwargs...)
	showfield(field, layer; ax=PyPlot.gca(), f=tiedrank, kwargs...)

Show geographic information and environmental information of a Mikrubi model.
The three principal components are reflexed in red, green, and blue. Keyword
argument `f = tiedrank` is a function acted on columns of `field.vars` as a 
whole. A possible alternative is `f = identity`.
"""
function showfield(field; ax=PyPlot.gca(), f=tiedrank, kwargs...)
	u(a) = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
	v(a) = u(f(a))
	r = field.dvar >= 1 ? v(field.vars[:, 1]) : fill(0.5, field.npixel)
	g = field.dvar >= 2 ? v(field.vars[:, 2]) : fill(0.5, field.npixel)
	b = field.dvar >= 3 ? v(field.vars[:, 3]) : fill(0.5, field.npixel)
	@assert size(field.locs, 2) >= 2
	x = field.locs[:, 1]
	y = field.locs[:, 2]
	ax.scatter(x, y, c=collect(zip(r, g, b)), s=2; kwargs...)
end
function showfield(field, layer; ax=PyPlot.gca(), f=tiedrank, kwargs...)
	u(a) = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
	v(a) = u(f(a))
	r = field.dvar >= 1 ? v(field.vars[:, 1]) : fill(0.5, field.npixel)
	g = field.dvar >= 2 ? v(field.vars[:, 2]) : fill(0.5, field.npixel)
	b = field.dvar >= 3 ? v(field.vars[:, 3]) : fill(0.5, field.npixel)
	xb, yb, _ = size(layer)
	matrix = fill(1., yb, xb, 3)
	for i = 1:field.npixel
		x, y = xy2ij(layer, field.locs[i, 1:2]...)
		matrix[y, x, :] .= r[i], g[i], b[i]
	end
	x, y = dims(layer)[1:2]
	x1, x2 = extrema(dimseq(x))
	y1, y2 = extrema(dimseq(y))
	handle = ax.imshow(matrix; 
		extent=(x1,x2, y1,y2), zorder=-1, aspect="auto", kwargs...)
	nna = minimum(matrix, dims=3) .< 1
	xval = any(nna, dims=1)[:]
	yval = any(nna, dims=2)[:]
	x3, x4 = minmax(dimseq(x)[findfirst(xval)], dimseq(x)[findlast(xval)+1])
	y3, y4 = minmax(dimseq(y)[findfirst(yval)], dimseq(y)[findlast(yval)+1])
	lim(z3, z4) = [1.05 -0.05; -0.05 1.05] * [z3, z4]
	ax.set_xlim(lim(x3, x4))
	ax.set_ylim(lim(y3, y4))
	return handle
end

"""
	showctpixels(ctpixels; ax=PyPlot.gca(), salt=20, kwargs...)
	showctpixels(ctpixels, layer; ax=PyPlot.gca(), salt=20, kwargs...)

Show a `Mikrubi.CtPixels`. Every county is assigned a hash color (influenced
by a fixed `salt` value also), and every pixel has the composite color from all
counties assigned to it. Empty cells are depicted white.
"""
showctpixels(ctpixels; ax=PyPlot.gca(), salt=20, kwargs...) = 
	showctpixels(ctpixels, ctpixels.indices; ax=ax, salt=salt, kwargs...)
function showctpixels(ctpixels, layer; ax=PyPlot.gca(), salt=20, kwargs...)
	xb, yb, _ = size(layer)
	ci = CartesianIndices((xb, yb))
	matrix = fill(1., yb, xb, 3)
	for (ct, pixel) = ctpixels.list
		x, y = ci[pixel].I
		matrix[y, x, :] .*= hashcolor(ct, salt=salt)
	end
	x, y = dims(layer)[1:2]
	x1, x2 = extrema(dimseq(x))
	y1, y2 = extrema(dimseq(y))
	handle = ax.imshow(matrix; 
		extent=(x1,x2, y1,y2), zorder=-1, aspect="auto", kwargs...)
	nna = minimum(matrix, dims=3) .< 1
	xval = any(nna, dims=1)[:]
	yval = any(nna, dims=2)[:]
	x3, x4 = minmax(dimseq(x)[findfirst(xval)], dimseq(x)[findlast(xval)+1])
	y3, y4 = minmax(dimseq(y)[findfirst(yval)], dimseq(y)[findlast(yval)+1])
	lim(z3, z4) = [1.05 -0.05; -0.05 1.05] * [z3, z4]
	ax.set_xlim(lim(x3, x4))
	ax.set_ylim(lim(y3, y4))
	return handle
end

function polygonline(geom)
	line = Tuple{Float64, Float64}[]
	parts, mat = geom2mat(geom)
	parts .+= 1
	for i = 1:size(mat,2)
		pt = (mat[:, i]...,)
		i in parts && push!(line, (NaN, NaN))
		push!(line, pt)
	end
	line
end

function shapelines(geoms)
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
	showshptable(shptable; ax=PyPlot.gca(), kwargs...)

Show lines from polygons in `shptable`. Identical segments are reduced as one.
"""
function showshptable(shptable; ax=PyPlot.gca(), kwargs...)
	shortline = shapelines(AG.getgeom.(shptable))
	ax.plot(first.(shortline), last.(shortline), "-k", lw=0.8; kwargs...)
end
