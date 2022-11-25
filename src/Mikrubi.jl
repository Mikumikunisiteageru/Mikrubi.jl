module Mikrubi

include("graphics.jl")

export Graphics

using GeoArrays

import ArchGDAL
import Printf: @sprintf
import DelimitedFiles: readdlm, writedlm
import Statistics: mean, std, cor
import MultivariateStats # fit, projection, PCA
import Optim: optimize, maximize, maximum, NelderMead, Newton, Options

export logistic, loglogistic
export readshape, goodcolumns, lookup, rasterize
export readlayers, writelayer, writelayers, makefield
export MikrubiField, readfield, writefield
export MikrubiModel, readmodel, writemodel
export readlist, writelist, fit
export predict, predictcounty, probcounties, samplecounties
export lipschitz

const MAXPCADIM = 4

# PROGRAMMING TRICKS AND GADGETS

"""
	I(a...)

Returns a tuple packaging all arguments, like a slurping version of `identity`.

# Example for broadcasting usage (somewhat resembling `zip`)
```julia
julia> I.(1, [2, 3, 4])
3-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (1, 3)
 (1, 4)

julia> I.([1, 0, -1], [2, 3, 4])
3-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (0, 3)
 (-1, 4)
```
"""
I(a...) = (a...,)

"""
	textwrap(str::AbstractString)

Gobbles all linefeeds ("`\\n`") inside `str` and replaces them with spaces
("` `"), so long strings can be wrapped to multiple lines in the codes, like 
the Python package `textwrap`.
"""
textwrap(str::AbstractString) = replace(str, r" *\n[ \t]*" => " ")

"""
	allsame(a::AbstractVector)

Returns `true` if all elements from `a` are identical, or otherwise `false`.

WARNING: An error is thrown if the vector `a` is empty.

# Examples
```julia
julia> allsame([1, 1, 2])
false

julia> allsame([1, 1, 1])
true

julia> allsame([1])
true

julia> allsame([])
ERROR: BoundsError: attempt to access 0-element Array{Any,1} at index [1]
Stacktrace:
 [1] getindex at .\\array.jl:787 [inlined]
 [2] allsame(::Array{Any,1}) at .\\REPL[9]:1
 [3] top-level scope at REPL[20]:1
```
"""
allsame(a::AbstractVector) = all(a[1] .== a[2:end])

"""
	colmatrix(vector::AbstractVector)
	colmatrix(matrix::AbstractMatrix)

Returns a one-column matrix if the argument is a vector, or the matrix itself
if the argument is already a matrix.
"""
colmatrix(vector::AbstractVector) = repeat(vector, 1, 1)
colmatrix(matrix::AbstractMatrix) = matrix

# MATHEMATIC TOOLS

"""
	logistic(x)

Computes `logistic(x) := 1 / (1 + e^x)`.
"""
logistic(x) = one(x) ./ (one(x) + exp(-x))

"""
	loglogistic(x)

Computes `log(logistic(x)) = -log(1 + e^x)`.
"""
loglogistic(x) = -log(one(x) + exp(-x))

# FUNCTIONS FOR RASTERIZING

"""
	interpolate(u1::Float64, u2::Float64, v1::Float64, v2::Float64)
	interpolate(u1, u2, v1, v2)

Returns an affine function as a closure, say `f`, which satisfies `f(u1) == v1`
and `f(u2) == v2`. 

WARNING: Arguments `u1` and `u2` cannot be identical. Beforehand judgement is
required. 
"""
interpolate(u1::Float64, u2::Float64, v1::Float64, v2::Float64) =
	u -> ((v1 - v2) * u + (u1*v2 - u2*v1)) / (u1 - u2)
interpolate(u1, u2, v1, v2) = interpolate(Float64.((u1, u2, v1, v2))...)

"""
	scanline(coords::AbstractMatrix{Float64}, parts::Vector{<:Integer}, xb, yb)

Performs scan line algorithm on polygon(s) defined by `coords` and `parts`, in
the rectangular field ranging from `0` to `xb` on longitude and from `0` to
`yb` on latitude, with grid resolution of `1` unit. Argument `coords` is a
matrix containing two rows representing longitude and latitude respectively,
and `parts` is a vector of integers denoting starting column index of each
ring in the polygon(s). When the center of a pixel lies within the polygon(s),
the center coordinates of the pixel are included in the output.
"""
function scanline(coords::AbstractMatrix{Float64}, 
		parts::AbstractVector{<:Integer}, xb, yb)
	buckets = Dict{Int, Vector{Float64}}()
	yc = round.(Int, coords[2, :])
	@inbounds for k = 1 : size(coords, 2)
		k in parts && continue
		yc[k-1] == yc[k] && continue
		y2x = interpolate(coords[2, k-1:k]..., coords[1, k-1:k]...)
		for y = min(yc[k-1], yc[k]) + 1 : max(yc[k-1], yc[k])
			push!(get!(buckets, y, Float64[]), y2x(y - 0.5))
		end
	end
	pixels = Tuple{Int,Int}[]
	@inbounds for y = keys(buckets)
		1 <= y <= yb || continue
		@assert iseven(length(buckets[y]))
		sort!(buckets[y])
		for k = 1:2:length(buckets[y])
			xs = max(1, round(Int, buckets[y][k]) + 1)
			xt = min(xb, round(Int, buckets[y][k + 1]))
			append!(pixels, I.(xs:xt, y))
		end
	end
	pixels
end

"""
	trimedge(x1, y1, x2, y2, xb, yb)

Trims the segment connecting `(x1, y1)` and `(x2, y2)` by a rectangle ranging
from `0` to `xb` on longitude and from `0` and `yb` on latitude.

# Example
```julia
julia> trimedge(-2.0, -3.0, 10.0, 6.0, 6.0, 7.0)
(2.0, 0.0, 6.0, 3.0)
```
"""
function trimedge(x1, y1, x2, y2, xb, yb)
	y1 > y2 && ((x1, y1, x2, y2) = (x2, y2, x1, y1))
	0 <= x1 <= xb && 0 <= x2 <= xb && 0 <= y1 && y2 <= yb &&
		return x1, y1, x2, y2
	(y2 < 0 || yb < y1 || max(x1, x2) < 0 || xb < min(x1, x2)) && 
		return NaN, NaN, NaN, NaN
	x1 == x2 && return x1, max(0, y1), x2, min(yb, y2)
	y1 == y2 && return max(0, min(x1, x2)), y1, min(xb, max(x1, x2)), y2
	x2y = interpolate(x1, x2, y1, y2)
	y2x = interpolate(y1, y2, x1, x2)
	yx = [ x2y(0) 0; x2y(xb) xb; 0 y2x(0); yb y2x(yb) ] # L, R, D, U
	yl, yh, xl, xh = sortslices(yx, dims=1)[2:3, :]
	yxf = [ y1 x1; yl xl; yh xh; y2 x2 ]
	yal, yah, xal, xah = sortslices(yxf, dims=1)[2:3, :]
	return xal, yal, xah, yah
end

"""
	strokeedge(x1, y1, x2, y2, xb, yb)

Gathers all unit pixels in the rectangle ranging from `0` to `xb` on longitude
and from `0` and `yb` on latitude that the segment connecting `(x1, y1)` and
`(x2, y2)` touches.
"""
function strokeedge(x1, y1, x2, y2, xb, yb)
	pixels = Tuple{Int,Int}[]
	x1, y1, x2, y2 = trimedge(x1, y1, x2, y2, xb, yb)
	isnan(x1) && return pixels
	xr = floor(Int, min(x1, x2)) + 1 : ceil(Int, max(x1, x2))
	yr = floor(Int, y1) + 1 : ceil(Int, y2)
	(length(xr) == 0 || length(yr) == 0) && return pixels
	(length(xr) == 1 || length(yr) == 1) && return I.(xr, yr)
	y2x = interpolate(y1, y2, x1, x2)
	if x1 < x2
		xleft = min(xb, max(0, x1))
		for y = yr[begin:end-1]
			xright = min(xb, max(0, y2x(y)))
			append!(pixels, I.(floor(Int, xleft) + 1 : ceil(Int, xright), y))
			xleft = xright
		end
		xright = min(xb, max(0, x2))
		return append!(pixels, 
			I.(floor(Int, xleft) + 1 : ceil(Int, xright), ceil(Int, y2)))
	else
		xright = min(xb, max(0, x1))
		for y = yr[begin:end-1]
			xleft = min(xb, max(0, y2x(y)))
			append!(pixels, I.(floor(Int, xleft) + 1 : ceil(Int, xright), y))
			xright = xleft
		end
		xleft = min(xb, max(0, x2))
		return append!(pixels, 
			I.(floor(Int, xleft) + 1 : ceil(Int, xright), ceil(Int, y2)))
	end
end

"""
	strokepath(coords::AbstractMatrix{Float64}, parts, xb, yb)

Gathers all unit pixels in the rectangle ranging from `0` to `xb` on longitude
and from `0` and `yb` on latitude that the path(s) defined by `coords` and
`parts` touch(es).
"""
function strokepath(coords::AbstractMatrix{Float64}, parts, xb, yb)
	x = coords[1, :]
	y = coords[2, :]
	pixels = Set{Tuple{Int,Int}}()
	@inbounds for k = 1 : size(coords, 2)
		k in parts && continue
		xr = range(max(1, floor(Int, min(x[k-1], x[k])) + 1),
		      stop=min(xb, ceil(Int, max(x[k-1], x[k]))))
		yr = range(max(1, floor(Int, min(y[k-1], y[k])) + 1),
		      stop=min(yb, ceil(Int, max(y[k-1], y[k]))))
		(length(xr) == 0 || length(yr) == 0) && continue
		if length(xr) == 1 || length(yr) == 1
			union!(pixels, I.(xr, yr))
			continue
		end
		union!(pixels, strokeedge(x[k-1], y[k-1], x[k], y[k], xb, yb))
	end
	pixels
end

"""
	geom2mat(polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon})
	geom2mat(multipolygon::ArchGDAL.IGeometry{ArchGDAL.wkbMultiPolygon})
	
Converts a `ArchGDAL.IGeometry` to a two-row matrix, with first row representing
the longitudes, and the second for latitudes.
"""
function geom2mat(polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon})
	parts = Int[]
	points = NTuple{3, Float64}[]
	for i = 0 : ArchGDAL.ngeom(polygon)-1
		ring = ArchGDAL.getgeom(polygon, i)
		push!(parts, length(points))
		for j = 0 : ArchGDAL.ngeom(ring)-1
			push!(points, ArchGDAL.getpoint(ring, j))
		end
	end
	n = length(points)
	mat = Matrix{Float64}(undef, 2, n)
	for i = 1:n
		mat[:, i] .= points[i][1:2]
	end
	return parts, mat
end
function geom2mat(multipolygon::ArchGDAL.IGeometry{ArchGDAL.wkbMultiPolygon})
	parts = Int[]
	points = NTuple{3, Float64}[]
	for k = 0 : ArchGDAL.ngeom(multipolygon)-1
		polygon = ArchGDAL.getgeom(multipolygon, k)
		for i = 0 : ArchGDAL.ngeom(polygon)-1
			ring = ArchGDAL.getgeom(polygon, i)
			push!(parts, length(points))
			for j = 0 : ArchGDAL.ngeom(ring)-1
				push!(points, ArchGDAL.getpoint(ring, j))
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

"""
	rasterize(geom::ArchGDAL.IGeometry, grid::GeoArray)
	
Returns coordinates of pixels of `grid` that `geom` the polygon(s) touch(es).
"""
function rasterize(geom::ArchGDAL.IGeometry, grid::GeoArray)
	@assert ArchGDAL.getgeomtype(geom) in 
		[ArchGDAL.wkbMultiPolygon, ArchGDAL.wkbPolygon]
	parts, mat = geom2mat(geom)
	parts = parts .+ 1
	coords = grid.f.linear \ (mat .- grid.f.translation)
	# yb, xb, _ = size(grid)
	xb, yb, _ = size(grid)
	pixels = scanline(coords, parts, xb, yb)
	union!(pixels, strokepath(coords, parts, xb, yb))
end

"""
	CtPixels

Alias for `Vector{Tuple{Tuple{Int, Int}, Int}}`, each element formatted as
`((i, j), ctid)`, where `i` and `j` are matrix index on row and column, and
`ctid` is identifier of a county.
"""
const CtPixels = Vector{Tuple{Tuple{Int, Int}, Int}}

getpixels(ctpixels::CtPixels) = first.(ctpixels)

getcounties(ctpixels::CtPixels) = last.(ctpixels)

"""
	rasterize(shptable::ArchGDAL.IFeatureLayer, grid::GeoArray)
	
Collects coordinates of pixels of `grid` that every polygon(s) from `shptable`
touch(es).
"""
function rasterize(shptable::ArchGDAL.IFeatureLayer, grid::GeoArray)
	geoms = ArchGDAL.getgeom.(shptable)
	ctpixels = CtPixels()
	@info textwrap("Rasterizing procedure may be slow. Please wait...")
	ndiv = 80
	unitlen = ndiv / length(geoms)
	print("-" ^ ndiv * "\r")
	for i = eachindex(geoms)
		print("#" ^ floor(Int, i * unitlen) * "\r")
		append!(ctpixels, I.(rasterize(geoms[i], grid), i))
	end
	println()
	ctpixels
end

# HELPER METHODS FOR SHAPEFILE HANDLING

"""
	readshape(path::AbstractString)

Reads a shape file (which ends with ".shp", ".geojson", or "gpkg") located at 
`path`.
"""
function readshape(path::AbstractString, index::Int=-1)
	ispath(path) || error(textwrap("The `path` given is invalid!"))
	if isdir(path)
		@info textwrap("The `path` is a directory. Now finding a random 
			file with extension `.shp`, `.geojson`, or `.gpkg` (from GADM) in 
			the directory by default.")
		files = readdir(path; join=true)
		isgeofile(f) = last(splitext(f)) in [".shp", ".geojson", ".gpkg"]
		i = findfirst(isgeofile, files)
		isnothing(i) && error(textwrap("No such a file in the directory!"))
		path = files[i]
	end
	dataset = ArchGDAL.read(path)
	if ArchGDAL.nlayer(dataset) == 1
		index = 0
	elseif ArchGDAL.nlayer(dataset) > 1 && index == -1
		display(dataset)
		error(textwrap("There are more than one data layers in the dataset! 
			Please designate the layer `index`."))
	end
	shptable = ArchGDAL.getlayer(dataset, index)
	isempty(shptable) && error(textwrap("The dataset does not contain any 
		shapes or geometries!"))
	return shptable
end

"""
	goodcolumns(shptable::ArchGDAL.IFeatureLayer)

Find all properties or fields of `shptable` where entries are all unique and 
either integers or strings (types where `isequal` is well-defined).
"""
function goodcolumns(shptable::ArchGDAL.IFeatureLayer)
	n = ArchGDAL.nfield(shptable)
	fields = Dict{String, Vector}()
	feature = first(shptable)
	for i = 0 : ArchGDAL.nfield(shptable)-1
		list = ArchGDAL.getfield.(shptable, i)
		if eltype(list) <: Union{Integer, AbstractString} && allunique(list)
			name = ArchGDAL.getname(ArchGDAL.getfielddefn(feature, i))
			fields[name] = list
		end
	end
	fields
end

@deprecate goodproperties(shptable) goodcolumns(shptable)

"""
	lookup(shptable::ArchGDAL.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entry)
	lookup(shptable::ArchGDAL.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entries::AbstractArray)
	lookup(shptable::ArchGDAL.IFeatureLayer)

Finds row(s) in the shape table whose `column` record(s) equal(s) to `entry` 
or elements of `entries`. When the third argument is an array, results are 
output as an array of the same shape by broadcasting. 
"""
function lookup(shptable::ArchGDAL.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entry)
	# if ! hasproperty(shptable, column)
	feature = first(shptable)
	i = ArchGDAL.findfieldindex(feature, column)
	if i == -1 || isnothing(i)
		columns = keys(goodcolumns(shptable))
		error(textwrap("No column in the shapefile named `$column`! 
			Recommended alternatives are $(join(repr.(columns), ", ")). 
			Please select one from them as `column`."))
	end
	shpcol = ArchGDAL.getfield.(shptable, i)
	index = findall(shpcol .== entry)
	if length(index) == 1
		return index[1]
	elseif length(index) > 1
		@warn textwrap("There are multiple entries in the column $column 
			equaling $(repr(entry))! All of these are returned, but please
			select exactly one from them for succeeding processing.")
		return index
	end # length(index) == 0 afterwards
	if isa(entry, AbstractString) && eltype(shpcol) <: Integer
		@warn textwrap("Types mismatched: `entry` is a string while the column
			designated contains integers. Please try again after using `parse`
			or `tryparse` to convert `entry` to `Int`. Nothing is returned
			here.")
		return nothing
	elseif isa(entry, Integer) && eltype(shpcol) <: AbstractString
		@warn textwrap("Types mismatched: `entry` is an integer while the
			column designated contains strings. Please try again after using 
			`string` to convert `entry`. Nothing is returned here.")
		return nothing
	end
	@warn textwrap("No matched record. Please check the input arguments 
		again. Nothing is returned here.")
	return nothing
end
lookup(shptable::ArchGDAL.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entries::AbstractArray) =
	lookup.([shptable], [column], entries)
lookup(shptable::ArchGDAL.IFeatureLayer) = lookup(shptable, "", 0)

# FUNCTIONS FOR HANDLING RASTER LAYERS

"""
	sortfilenames!(filenames::AbstractVector{<:AbstractString})

Identifies the distinctive parts among `filenames`, and sorts `filenames`
according to the order of those parts. If all of the distinctive parts are
decimal numerals, they are sorted as integers.

# Examples
```julia
julia> sortfilenames!(["bio_9.tif", "bio_10.tif", "bio_1.tif"])
[ Info: Totally 3 raster files recognized in the directory, in the form of bio_*.tif, where the asterisk means 1, 9, 10.
3-element Array{String,1}:
 "bio_1.tif"
 "bio_9.tif"
 "bio_10.tif"

julia> sortfilenames!(["bio_09.tif", "bio_10.tif", "bio_01.tif"])
[ Info: Totally 3 raster files recognized in the directory, in the form of bio_*.tif, where the asterisk means 01, 09, 10.
3-element Array{String,1}:
 "bio_01.tif"
 "bio_09.tif"
 "bio_10.tif"
```
"""
function sortfilenames!(filenames::AbstractVector{<:AbstractString})
	if length(filenames) == 0
		error(textwrap("No available raster file(s) in the directory!"))
	elseif length(filenames) == 1
		@info textwrap("Only one available raster file recognized in the 
			directory, namely $(splitpath(filenames[1])[end]).")
		return
	end
	s = 1
	while allsame([fn[s] for fn = filenames])
		s += 1
	end
	t = 0
	while allsame([fn[end-t] for fn = filenames])
		t += 1
	end
	fileids = [fn[s:end-t] for fn = filenames]
	perm = sortperm(fileids)
	try
		perm = sortperm(parse.(Int, fileids))
	catch ; end
	@info textwrap("Totally $(length(filenames)) raster files recognized in
		the directory, in the form of 
		$(splitpath(filenames[1][1:s-1])[end])*$(filenames[1][end-t+1:end]),
		where the asterisk means $(join(fileids[perm], ", ")).")
	filenames .= filenames[perm]
end

"""
	copylayer(layer::GeoArray; 
		fixmissing=true, verylargebase=1_000_000_000_000)
	
Copies a `GeoArray` layer (so that the matrix is not read-only).

By setting `fixmissing = true`, all pixels with absolute values greater than
`verylargebase` (= 10^12 by default) are identified as `missing`.
"""
function copylayer(layer::GeoArray; 
		fixmissing=true, verylargebase=1_000_000_000_000)
	if fixmissing
		A = Array{Union{Missing, eltype(layer.A)}}(layer.A[:, :, :])
		idx = .!ismissing.(A)
		absA = abs.(A[idx])
		verylarge = max(verylargebase, maximum(absA))
		A[findall(idx)[absA .>= verylarge]] .= missing
	else
		A = copy(layer.A)
	end
	GeoArray(A, layer.f, layer.crs)
end

"""
	readlayers(dir::AbstractString; extset=nothing)

Returns all raster layers as a vector of `GeoArray`s from the directory `dir`.
All elements in the vector returned are modifiable, and modifications do
not affect the files on the disk. 

By setting `extset` to a container of extensions (e.g. `Set(".tif")`, or
`[".tiff"]`), only files ending with its elements are read. Default value
is `nothing`, which means no file extension filtering.
"""
function readlayers(dir::AbstractString; extset=nothing)
	if isnothing(extset)
		check = isfile
	else
		check = path -> isfile(path) && last(splitext(path)) in extset
	end
	filenames = filter(check, readdir(dir, join=true))
	sortfilenames!(filenames)
	layers = GeoArrays.read.(filenames)
	for i = 2:length(layers)
		(size(layers[i]) != size(layers[1]) || layers[1].f != layers[2].f) &&
			error(textwrap("Resolutions, shapes, or ranges of the layers are
				not consistent!"))
	end
	copylayer.(layers, fixmissing=true)
end

"""
	writelayer(path::AbstractString, layer::GeoArray)

Writes `layer` to the disk at `path`. Alias for `GeoArrays.write!`.
"""
writelayer(path::AbstractString, layer::GeoArray) =
	GeoArrays.write!(path, layer)

"""
	writelayers(paths::AbstractVector{<:AbstractString}, 
		layers::AbstractVector)
	writelayers(pathformula::AbstractString, layers::AbstractVector)

Writes `layers` to `paths` respondingly, or a series of paths generated by the
`pathformula` where an asterisk is used for wildcard and replaced by numbers.
"""
function writelayers(paths::AbstractVector{<:AbstractString}, 
		layers::AbstractVector)
	length(paths) != length(layers) &&
		error(textwrap("Arguments `paths` and `layers` must have the same 
			length as vectors!"))
	for (path, layer) = zip(paths, layers)
		writelayer(path, layer)
	end
end
function writelayers(pathformula::AbstractString, layers::AbstractVector)
	occursin("*", pathformula) ||
		error(textwrap("Argument `pathformula` must contain an asterisk 
			(\"`*`\") as wildcard!"))
	for i = 1:length(layers)
		path = replace(pathformula, "*" => i)
		writelayer(path, layers[i])
	end
end

"""
	masklayers!(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels)

Masks the `layers` in a way that only pixels present in `ctpixels` are kepts,
while all other unmentioned pixels are set to `missing`.
"""
function masklayers!(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels)
	xb, yb, nband = size(layers[1])
	mask = trues(xb, yb, nband)
	for (x, y) = getpixels(ctpixels)
		mask[x, y, :] .= false
	end
	for i = eachindex(layers)
		layers[i].A[mask] .= missing
	end
	layers
end

"""
	extractlayers(layers::AbstractVector{<:GeoArray})

Extracts the non-`missing` pixels from the `layers`, and combines them into a
matrix, whose rows representing pixels and columns representing variables.

`extractlayers` is the inverse function of `makelayers`.
"""
function extractlayers(layers::AbstractVector{<:GeoArray})
	nband = last.(size.(layers))
	all(nband .== 1) || length(nband) == 1 ||
		@error textwrap("Multiple raster objects containing multiple bands! 
			Please reorganize the raster layers either in only one file as 
			multiple bands or in multiple files each of which containing only 
			one band.")
	if all(nband .== 1)
		variables = [layer.A[:, :, 1] for layer = layers]
	else # length(nband) == 1
		variables = [layers[1].A[:, :, i] for i = nband[1]]
	end
	idx = ismissing.(variables[1])
	missingcount = sum(idx)
	for i = 2:length(variables)
		idx .|= ismissing.(variables[i])
	end
	missingcount != sum(idx) && 
		@warn textwrap("Missing-data pixels are inconsistent across the 
			layers! Only pixels of validate values across all raster layers are
			considered in calculation.")
	hcat([var[.!idx] for var = variables]...), findall(.!idx)
end

"""
	emptylayer(grid::GeoArray)

Creates a `GeoArray` with the same size and other attributes as `grid`, but
full of `missing`s in its matrix.
"""
function emptylayer(grid::GeoArray)
	xb, yb, _ = size(grid)
	A = Array{Union{Missing, eltype(grid.A)}}(fill(missing, xb, yb, 1))
	GeoArray(A, grid.f, grid.crs)
end

"""
	makelayers(matrix::AbstractMatrix, idx::AbstractVector, grid::GeoArray)

Returns a vector of layers according to the grid from `grid` and content values
from columns of `matrix`. See also `makelayer`.

`makelayers` is the inverse function of `extractlayers`.
"""
function makelayers(matrix::AbstractMatrix, 
		idx::AbstractVector, grid::GeoArray)
	npixel, mvar = size(matrix)
	npixel == length(idx) ||
		error(textwrap("The number of columns in `matrix` has to be identical
			to the length of `idx`!"))
	layers = [emptylayer(grid) for _ = 1:mvar]
	for i = 1:mvar
		layers[i].A[idx] .= matrix[:, i]
	end
	layers
end

"""
	makelayer(vector::AbstractVector, idx::AbstractVector, grid::GeoArray)

Returns a layer according to the grid from `grid` and content values from 
`vector`. See also `makelayers`.
"""
makelayer(vector::AbstractVector, idx::AbstractVector, grid::GeoArray) =
	makelayers(colmatrix(vector), idx, grid)[1]

"""
	dftraverse!(beststate, bestscore, state, score, depth, maxdepth, 
		incompat, scoremat)

Finds the index combination that
- firstly containing as many indices as possible, and
- secondly with the lowest pairwise sum from submatrix of `scoremat`,
such that no indices `i` and `j` coexist as long as `incompat[i][j] == true`.

The result is stored as the only element of `beststate`, with its score
decided by the two criteria above stored as the only element of `bestscore`.

# Example
```julia
julia> beststate = Vector(undef, 1);

julia> bestscore = [(0, 0.0)];

julia> dftraverse!(beststate, bestscore, Int[], (0, 0.0), 1, 3,
           Bool[0 0 1; 0 0 0; 1 0 0],
           [0.0 0.6 0.3; 0.6 0.0 0.9; 0.3 0.9 0.0]);

julia> beststate
1-element Array{Any,1}:
 [1, 2]

julia> bestscore
1-element Array{Tuple{Int64,Float64},1}:
 (2, -0.6)
```
"""
function dftraverse!(beststate, bestscore, state, score, depth, maxdepth, 
		incompat, scoremat)
	if depth > maxdepth	
		if score > bestscore[1]
			bestscore[1] = score
			beststate[1] = state
		end
		return
	end
	if ! any(incompat[state, depth])
		dftraverse!(beststate, bestscore, vcat(state, depth), 
			score .+ (1, -sum(scoremat[state, depth])), depth + 1, maxdepth,
			incompat, scoremat)
	end
	dftraverse!(beststate, bestscore, state, score, depth + 1, maxdepth,
		incompat, scoremat)
end

"""
	selectvars(matrix::Matrix, rabsthres=0.8)

Selects as many variables as possible from `matrix` such that no pairwise
Pearson coefficient among them exceeds `rabsthres` and their sum is minimal.

# Example
```julia
julia> selectvars([1. 4. 7.; 2. 5. 8.; 3. 9. 27.], rabsthres=0.9)
2-element Array{Int64,1}:
 1
 3
```
"""
function selectvars(matrix::Matrix; rabsthres=0.8)
	rabsmat = abs.(cor(matrix, dims=1))
	incompat = rabsmat .> rabsthres
	beststate = Vector(undef, 1)
	bestscore = [(0, zero(eltype(matrix)))]
	dftraverse!(beststate, bestscore, Int[], bestscore[1], 1, size(matrix, 2),
		incompat, rabsmat)
	beststate[1]
end

"""
	princompvars(submat::Matrix; nprincomp=3)

Performs principal component analysis on `submat` whose columns represents
variables, and combines the `nprincomp` principal components into a matrix, and
returns the result matrix as well as the affine transformation
`(colmean, projwstd)`, such that
the result matrix == `(submat .- colmean) * projwstd`.
"""
function princompvars(submat::Matrix; nprincomp=3)
	nprincomp > MAXPCADIM && 
		@warn textwrap("It is strongly recommended that no more than four
			principal components are used for Mikrubi, or parameter space would
			be highly ill-conditioned!")
	colmean = mean(submat, dims=1)
	colstd = std(submat, dims=1)
	minimum(colstd) .<= 1e-11 &&
		error(textwrap("Some variable among the layers has a very small
			deviation (in other words, is (nearly) constant), which directly
			leads to ill condition of succeeding calculations."))
	submat01 = (submat .- colmean) ./ colstd
	pca = MultivariateStats.fit(MultivariateStats.PCA, 
		submat01', maxoutdim=nprincomp)
	pcadim = size(pca)[2]
	pcadim < nprincomp &&
		@warn textwrap("Only $(size(pca)[2]) principal component(s) (less than
			`nprincomp` = $nprincomp) is/are used to expressed the selected
			layers in the principal component analysis.")
	proj = MultivariateStats.projection(pca)
	projwstd = proj ./ colstd[:]
	submat01 * proj, (colmean, projwstd)
end

"""
	makefield(layers::AbstractVector{<:GeoArray}, 
		shptable::ArchGDAL.IFeatureLayer; rabsthres=0.8, nprincomp=3)
	makefield(layers::AbstractVector{<:GeoArray}, 
		shptable::ArchGDAL.IFeatureLayer, players::Vector{<:GeoArray}; 
		rabsthres=0.8, nprincomp=3)
	makefield(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels; 
		rabsthres=0.8, nprincomp=3)
	makefield(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels,
		players::Vector{<:GeoArray}; rabsthres=0.8, nprincomp=3)

Generates a Mikrubi field as well as processed variable layers from `layers` 
and `shptable` or `ctpixels`, by 
0. (rasterizing the `shptable` to `ctpixels` using `rasterize`,)
1. masking the `layers` with `ctpixels` (using `Mikrubi.masklayers!`),
2. extracting non-missing pixels from `layers` (using `Mikrubi.extractlayers`),
3. selecting less correlated variables (using `Mikrubi.selectvars`), and
4. doing the principal component analysis (using `Mikrubi.princompvars`).

# Optional keyword arguments

- `rabsthres`: threshold of collinearity.
Absolute value of Pearson correlation efficient greater than this threshold 
is identified as collinearity and the two variables are thus incompatible. 
- `nprincomp`: expected number of principal components of the variables.

# Notes about `players`

when `players` is present in the argument list, raster layers it contains
experience the same process including subsetting, selecting, and taking
principal components, and results are packed and returned in the third place.
Users must assure that `players` has the same length as `layers`, and
their elements are corresponding in order. This would be useful when the
prediction is in another geographic range or at another time.
"""
makefield(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels; 
	rabsthres=0.8, nprincomp=3) = 
		_makefield(layers, ctpixels, rabsthres, nprincomp)[1:2]
function makefield(layers::AbstractVector{<:GeoArray}, ctpixels::CtPixels,
		players::Vector{<:GeoArray}; rabsthres=0.8, nprincomp=3)
	length(layers) == length(players) ||
		error(textwrap("Vectors `players` must have the same length as 
			`layers`!"))
	field, elayers, subindex, colmean, projwstd =
		_makefield(layers, ctpixels, rabsthres, nprincomp)
	pmatrix, pidx = extractlayers(players)
	subpmat = pmatrix[:, subindex]
	epmat = (subpmat .- colmean) * projwstd
	eplayers = makelayers(epmat, pidx, players[1])
	return field, elayers, eplayers
end
function _makefield(layers, ctpixels, rabsthres, nprincomp)
	layers = copylayer.(layers)
	masklayers!(layers, ctpixels)
	matrix, idx = extractlayers(layers)
	subindex = selectvars(matrix)
	submat = matrix[:, subindex]
	emat, (colmean, projwstd) = princompvars(submat, nprincomp=nprincomp)
	field = MikrubiField(ctpixels, idx, emat, layers[1])
	elayers = makelayers(emat, idx, layers[1])
	field, elayers, subindex, colmean, projwstd
end
makefield(layers::AbstractVector{<:GeoArray}, shptable::ArchGDAL.IFeatureLayer;
	rabsthres=0.8, nprincomp=3) =
		makefield(layers, rasterize(shptable, layers[1]),
			rabsthres=rabsthres, nprincomp=nprincomp)
makefield(layers::AbstractVector{<:GeoArray}, shptable::ArchGDAL.IFeatureLayer,
	players::AbstractVector{<:GeoArray}; rabsthres=0.8, nprincomp=3) =
		makefield(layers, rasterize(shptable, layers[1]),
			players, rabsthres=rabsthres, nprincomp=nprincomp)

# THE CORE OF MIKRUBI

"""
	dvar2dparam(dvar::Int)

Transforms dimensionality of an environmental space to the dimensionality of
the induced parameter space, i.e., returns the degrees of freedom for 
positive-definite quadratic functions mapping a `dvar`-dimensional linear space
into real numbers.

# Examples
```julia
julia> dvar2dparam(1)
3

julia> dvar2dparam(3)
10
```
"""
dvar2dparam(dvar::Int) = ((dvar + 1) * (dvar + 2)) >> 1

"""
	MikrubiField(ctids, locs, vars)

Constructs a Mikrubi field containing a number of pixels or points, using the
arguments which should always have the same number of rows
- `ctids::Vector`: a vector containing the county identifiers
- `locs::Array{<:Real}`: an array of geographic coordinates
- `vars::Matrix{<:AbstractFloat}`: an array of environmental coordinates
"""
struct MikrubiField{T, U <: Real, V <: AbstractFloat}
	ctids::Vector{T}
	locs::VecOrMat{U}
	vars::Matrix{V}
	npixel::Int
	mcounty::Int
	ids::Vector{T}
	starts::Dict{T, Int}
	stops::Dict{T, Int}
	dvar::Int
	function MikrubiField(ctids::AbstractVector{T}, locs::AbstractVecOrMat{U}, 
			vars::AbstractMatrix{V}) where {T, U <: Real, V <: AbstractFloat}
		ctids, locs, vars = copy(ctids), copy(locs), copy(vars)
		npixel = length(ctids)
		npixel == size(locs, 1) == size(vars, 1) ||
			error(textwrap("Arguments `ctids`, `locs`, and `vars`
				should always have the same number of rows!"))
		if ! issorted(ctids)
			perm = sortperm(ctids)
			ctids, locs, vars = ctids[perm], locs[perm, :], vars[perm, :]
		end
		ids = unique(ctids)
		mcounty = length(ids)
		size(vars, 2) >= 5 && 
			@warn textwrap("It is strongly recommended that no more than four
				principal components are used for Mikrubi, or parameter space
				would be highly ill-conditioned!")
		starts = Dict(reverse(ctids) .=> npixel:-1:1)
		stops  = Dict(        ctids  .=> 1:npixel   )
		dvar = size(vars, 2)
		new{T, U, V}(ctids, locs, vars, 
			npixel, mcounty, ids, starts, stops, dvar)
	end
end
MikrubiField(ctids, locs, vars) =
	MikrubiField(ctids, colmatrix(locs), colmatrix(float.(vars)))

"""
	MikrubiField(ctpixels::CtPixels, idx::Vector, 
		projmat::Matrix, grid::GeoArray)

Constructs a `MikrubiField` from `ctpixels`, `idx`, `projmat`, and `grid`. See also 
`makefield`.
"""
function MikrubiField(ctpixels::CtPixels, idx::AbstractVector, 
		projmat::AbstractMatrix, grid::GeoArray)
	ctids = getcounties(ctpixels)
	npixel = length(ctids)
	mvar = size(projmat, 2)
	revidx = Dict(getfield.(idx, :I) .=> 1:length(idx))
	locs = Matrix{Float64}(undef, npixel, 2)
	vars = Matrix{eltype(projmat)}(undef, npixel, mvar)
	found = trues(npixel)
	for i = 1:npixel
		coord = first(ctpixels[i])
		if haskey(revidx, coord)
			locs[i, :] .= centercoords(grid, collect(coord))
			vars[i, :] .= projmat[revidx[coord], :]
		else
			found[i] = false
		end
	end
	sumfound = sum(found)
	ratiofound = sumfound / npixel
	0.9 < ratiofound < 1 &&
		@info textwrap("Among the $npixel pixels, $(npixel-sumfound) 
			($(@sprintf("%.1f", 100 * (1 - ratiofound)))%) is/are discarded 
			for lacking values.")
	ratiofound <= 0.9 &&
		@warn textwrap("Among the $npixel pixels, $(npixel-sumfound) 
			($(@sprintf("%.1f", 100 * (1 - ratiofound)))%) is/are discarded 
			for lacking values!")
	MikrubiField(ctids[found], locs[found, :], vars[found, :])
end

function Base.show(io::IO, field::MikrubiField)
	print(io, textwrap("Mikrubi Field: 
		geo_dim = $(size(field.locs, 2)),
		env_dim = $(field.dvar),
		$(field.npixel) pixels,
		and $(field.mcounty) counties"))
end

"""
	writefield(path::AbstractString, field::MikrubiField)

Writes a Mikrubi field to file.
"""
function writefield(path::AbstractString, field::MikrubiField)
	headerstring = "I" * "L" ^ size(field.locs, 2) * "V" ^ field.dvar
	header = string.(hcat(headerstring...))
	body = hcat(Array{Any}(field.ctids), field.locs, field.vars)
	writedlm(path, vcat(header, body))
end

"""
	readfield(path::AbstractString)

Reads a Mikrubi field from file.
"""
function readfield(path::AbstractString)
	body, header = readdlm(path, Any, header=true)
	heads = header[:]
	Set(heads) == Set(["I", "L", "V"]) && findlast(heads .== "I") == 1 &&
		findfirst(heads .== "V") - findlast(heads .== "L") == 1 ||
			error(textwrap("The file at `path` is not a well-formatted file!"))
	ctids = [body[:, heads .== "I"]...]
	locs = Real.(body[:, heads .== "L"])
	vars = AbstractFloat.(body[:, heads .== "V"])
	MikrubiField(ctids, locs, vars)
end

"""
	MikrubiModel(dvar::Int, params::Vector{<:AbstractFloat})

Constructs a Mikrubi Model from a dimensionality `dvar` and a parameter vector
`params`. The equation must hold for `dvar2dparam(dvar) == length(params)`.

Mikrubi Models can be obtained from the function `fit`, and can be used in the
function `predict`.
"""
struct MikrubiModel{V <: AbstractFloat}
	dvar::Int
	params::Vector{V}
	function MikrubiModel(dvar::Int, 
			params::Vector{V}) where {V <: AbstractFloat}
		dvar2dparam(dvar) != length(params) &&
			error(textwrap("Length of `params` and value of `dvar` are
				incompatible! The following equation must hold:
				`(dvar+1) * (dvar+2)) / 2 == length(params)`."))
		new{V}(dvar, params)
	end
end

"""
	writemodel(path::AbstractString, model::MikrubiModel)

Writes a Mikrubi model to file.
"""
writemodel(path::AbstractString, model::MikrubiModel) =
	writedlm(path, vcat(Any[model.dvar], model.params))

"""
	readmodel(path::AbstractString)

Reads a Mikrubi model from file.
"""
function readmodel(path::AbstractString)
	vector = readdlm(path, header=false)[:]
	isinteger(vector[1]) &&
		length(vector)-1 == dvar2dparam(Int(vector[1])) ||
			error(textwrap("The file at `path` is not a well-formatted file!"))
	MikrubiModel(Int(vector[1]), vector[begin+1:end])
end

"""
	writemodel(path::AbstractString, list::AbstractVector)

Writes any list or vector to file.
"""
writelist(path::AbstractString, list::AbstractVector) = writedlm(path, list)

"""
	readlist(path::AbstractString)

Reads any list of vector from file.
"""
readlist(path::AbstractString) = [readdlm(path, Any, header=false)...]

"""
	decomparams(p::AbstractVector, d::Int)
	decomparams(model::MikrubiModel)

Returns parameter decomposition `At`, `b`, `c`, where 
- `At` is a lower triangular matrix of size `(d, d)`,
- `b` is a column vector of size `d`, and
- `c` is a scalar.

WARNING: The vector `p` must have length `dvar2dparam(d)`.

# Example
```julia
julia> decomparams(collect(1:10), 3)
([1 0 0; 2 3 0; 4 5 6], [7, 8, 9], 10)
```
"""
function decomparams(p::AbstractVector, d::Int)
	At = zeros(eltype(p), d, d)
	i = j = k = 1
	while i <= d
		At[i, j] = p[k]
		i == j ? i += j = 1 : j += 1
		k += 1
	end
	At, p[k:k+d-1], p[k+d]
end
decomparams(model::MikrubiModel) = decomparams(model.params, model.dvar)

"""
	loglike(field::MikrubiField, params::AbstractVector)

Computes the log-likelihood of being absent in each pixel given `field` and
`params`.
"""
function loglike(field::MikrubiField, params::AbstractVector)
	At, b, c = decomparams(params, field.dvar)
	pv = field.vars
	loglogistic.(sum((pv * At) .^ 2, dims=2)[:] .+ pv * b .+ c)
end

"""
	findnearest(loc::AbstractVecOrMat{<:Real}, field::MikrubiField)

Returns the row number in `field.locs` which is the nearest to the given 
coordinates.
"""
function findnearest(loc::AbstractVecOrMat{<:Real}, field::MikrubiField)
	loc = loc[:]'
	length(loc) == size(field.locs, 2) ||
		error(textwrap("Dimensionality of `loc` ($(length(loc))) is 
			incompatible with the geographic dimensionality of `field` 
			($(size(field.locs, 2)))!"))
	eucldist2 = sum((loc .- field.locs) .^ 2, dims=2)[:]
	findmin(eucldist2)[2]
end

"""
	findnearests(loc::AbstractVector{<:AbstractVecOrMat}, field::MikrubiField)
	findnearests(loc::AbstractMatrix{<:Real}, field::MikrubiField)

Returns the row numbers in `field.locs` which are the nearest to each of the 
given coordinates. Duplicate results are reduced to one.
"""
findnearests(loc::AbstractVector{<:AbstractVecOrMat}, field::MikrubiField) =
	unique(findnearest.(loc, [field]))
findnearests(loc::AbstractMatrix, field::MikrubiField) =
	unique(findnearest.(eachrow(float.(loc)), [field]))

"""
	energy(field::MikrubiField, counties, params::AbstractVector)
	energy(vars::AbstractMatrix, params::AbstractVector)

Computes the opposite log-likelihood that the occupied counties or occupied 
coordinates are sampled.
The result is taken opposite sign for optimization, and therefore the function
is called `energy`.
"""
function energy(field::MikrubiField, counties, params::AbstractVector)
	e = loglike(field, params)
	for o = counties
		start = field.starts[o]
		stop = field.stops[o]
		subsum = sum(e[start:stop])
		e[start:stop] .= 0
		e[start] = log(1 - exp(subsum))
	end
	- sum(e)
end
function energy(vars::AbstractMatrix, params::AbstractVector)
	At, b, c = decomparams(params, size(vars, 2))
	prallp = 1 .- logistic.(sum((vars * At) .^ 2, dims=2)[:] .+ vars * b .+ c)
	- sum(log.(prallp))
end

"""
	fit(field::MikrubiField, counties, coords=zeros(0, 0); 
		optresult=[], iterations=3_000_000, kwargs...)

Numerically finds the Mikrubi model maximizing the likelihood that the occupied
counties as well as the occupied coordinates are sampled in the given Mikrubi 
field. The optimization result is stored in the container `optresult` for 
debugging.
"""
function fit(field::MikrubiField, counties, coords=zeros(0, 0); 
		optresult=[], iterations=39000, kwargs...)
	valcounties = intersect(counties, field.ids)
	indcoords = findnearests(coords, field)
	isempty(valcounties) && isempty(indcoords) &&
		error(textwrap("No meaningful occupied counties or coordinates!"))
	cdvars = field.vars[indcoords, :]
	fun(params) = energy(field, valcounties, params) + energy(cdvars, params)
	# zeroes = zeros(eltype(field.vars), dvar2dparam(field.dvar))
	# TODO: fix the bug in Optim.jl involving NaN32
	zeroes = zeros(Float64, dvar2dparam(field.dvar))
	@info textwrap("Now minimizing the opposite likelihood function...")
	result = optimize(fun, zeroes, NelderMead(), Options(
		iterations=iterations, show_trace=true, show_every=500; kwargs...))
	result.iteration_converged && 
		@warn textwrap("The optimizing routine has reached the maximum 
			iteration count (`iterations = $iterations`), and thus the 
			maximizer may be unreliable. Please try to enlarge the parameter.")
	push!(optresult, result)
	@info textwrap("Maximized log-likeliness: $(-result.minimum)")
	# MikrubiModel(field.dvar, result.minimizer)
	MikrubiModel(field.dvar, eltype(field.vars).(result.minimizer))
end

"""
	predict(matrix::AbstractMatrix, model::MikrubiModel)
	predict(layers::AbstractVector{<:GeoArray}, model::MikrubiModel)
	predict(field::MikrubiField, model::MikrubiModel)

Predicts the probability of presence according to processed climatic factors
(`matrix` / `layers`) or on the Mikrubi field.
"""
function predict(matrix::AbstractMatrix, model::MikrubiModel)
	size(matrix, 2) == model.dvar ||
		error(textwrap("The number of columns of the matrix
			($(size(matrix, 2))) is different from the dimensionality of the
			model ($(model.dvar))!"))
	At, b, c = decomparams(model)
	1 .- logistic.(sum((matrix * At) .^ 2, dims=2)[:] .+ matrix * b .+ c)
end
function predict(layers::AbstractVector{<:GeoArray}, model::MikrubiModel)
	matrix, idx = extractlayers(layers)
	makelayer(predict(matrix, model), idx, layers[1])
end
predict(field::MikrubiField, model::MikrubiModel) =
	Dict(id => predictcounty(field, model, id) for id = field.ids)

"""
	predictcounty(field::MikrubiField, model::MikrubiModel, county)

Returns the geographic coordinates of pixels in the county sorted by the
likeliness of being occupied.
"""
function predictcounty(field::MikrubiField, model::MikrubiModel, county)
	county in field.ids ||
		error(textwrap("The Mikrubi field has no such county!"))
	idx = field.starts[county] : field.stops[county]
	probs = predict(field.vars[idx, :], model)
	locs = collect.(eachrow(field.locs[idx, :]))
	sort!(I.(locs, probs), rev=true, by=last)
end

"""
	probpixels(field::MikrubiField, model::MikrubiModel)

Returns the probability of being occupied for every pixel in the `field`. 
"""
probpixels(field::MikrubiField, model::MikrubiModel) = 
	1 .- exp.(loglike(field, model.params))

"""
	probcounties(field::MikrubiField{T, U, V}, model::MikrubiModel{V})

Returns the probability of being occupied for every county in the `field`. 
"""
function probcounties(field::MikrubiField{T, U, V}, 
		model::MikrubiModel{V}) where {T, U <: Real, V <: AbstractFloat}
	ll = loglike(field, model.params)
	logmprcounties = Dict{T, V}()
	for i = 1:field.npixel
		id = field.ctids[i]
		logmprcounties[id] = ll[i] + get(logmprcounties, id, 0)
	end
	prcounties = Dict{T, V}()
	for id = field.ids
		prcounties[id] = 1 - exp(logmprcounties[id])
	end
	prcounties
end

"""
	samplecounties(field::MikrubiField, model::MikrubiModel)

Samples the counties according to their probability of being occupied.
"""
function samplecounties(field::MikrubiField, model::MikrubiModel)
	prcounties = probcounties(field, model)
	@inline bernoulli(p) = rand() <= p
	sort!([k for k = keys(prcounties) if bernoulli(prcounties[k])])
end

"""
	loglipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false)

Calculates the (logarithmic) maximum gradient (in norm) of the probability of 
presence over the `field`. When `wholespace=false` (default), the maximum is 
taken among the points contained in `field`; otherwise it is taken around the 
whole space.
"""
function loglipschitz(model::MikrubiModel, 
		field::MikrubiField; wholespace=false)
	At, b, c = decomparams(model)
	biM = 2 * At * At'
	q(z) = sum((z*At) .^ 2) + sum(z * b) + c
	llpll(qz) = loglogistic(qz) + loglogistic(-qz)
	loglip(z) = llpll(q(z')) + log(sum((biM * z + b).^2))/2
	vars = field.vars
	m, id = findmax(loglip.(eachrow(vars)))
	if wholespace == false
		return m
	else
		return maximum(maximize(loglip, vars[id, :], 
			model.dvar == 1 ? Newton() : NelderMead()))
	end
end

"""
	lipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false)

Calculates the maximum gradient (in norm) of the probability of presence over 
the `field`. When `wholespace=false` (default), the maximum is taken among the 
points contained in `field`; otherwise it is taken around the whole space.
"""
lipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false) = 
	exp(loglipschitz(model, field; wholespace=wholespace))

end # module
