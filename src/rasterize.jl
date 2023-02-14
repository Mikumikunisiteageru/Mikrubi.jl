# src/prepare.jl

"""
	CtPixels
	
	CtPixels(indices::Raster{Int})

Collector for county-specific rasterization results, whose `list` contains 
county-pixel tuples. Can only be instantiated from an index raster 
(see [`indicate`](@ref)).
"""
struct CtPixels
	indices::Raster{Int}
	list::Vector{Tuple{Int,Int}}
	CtPixels(indices::Raster{Int}) = new(indices, Tuple{Int,Int}[])
end

Base.length(ctpixels::CtPixels) = length(ctpixels.list)

"""
	getpixels(ctpixels::CtPixels) :: Vector{Int}

Get pixel indices from `ctpixels`.
"""
getpixels(ctpixels::CtPixels) = last.(ctpixels.list)

"""
	getcounties(ctpixels::CtPixels) :: Vector{Int}

Get county indices from `ctpixels`.
"""
getcounties(ctpixels::CtPixels) = first.(ctpixels.list)

"""
	getpixel(ctpixels::CtPixels, i::Int) :: Int

Get the pixel index of the `i`-th county-pixel tuple in `ctpixels`.
"""
getpixel(ctpixels::CtPixels, i::Int) = last(ctpixels.list[i])

"""
	getcounty(ctpixels::CtPixels, i::Int) :: Int

Get the county index of the `i`-th county-pixel tuple in `ctpixels`.
"""
getcounty(ctpixels::CtPixels, i::Int) = first(ctpixels.list[i])

"""
	indicate(layer::Raster) :: Raster{Int}

Build an index raster `indices` from `layer`. The value of an array element in 
`indices` is either (1) `0` for missing, if the corresponding element in 
`layer` is missing; or (2) the integer index of the array element, otherwise. 
"""
function indicate(layer::Raster)
	indices = Int.(zero(layer))
	indices[:] .= eachindex(indices)
	indices[.!boolmask(layer)] .= 0
	return rebuild(indices, missingval=0)
end

"""
	register!(ctpixels::CtPixels, ct::Int, pixel::Int) :: Int

Push a county-pixel tuple into `ctpixels`, if `pixel` is not zero. For 
convenience, the value of `pixel` is returned. 
"""
function register!(ctpixels::CtPixels, ct::Int, pixel::Int)
	iszero(pixel) || push!(ctpixels.list, (ct, pixel))
	return pixel
end

"""
	register!(ctpixels::CtPixels, ct::Int) :: Function

Create a function that accepts a `pixel`, pushes the county-pixel tuple into 
`ctpixels`, and finally returns the value of `pixel`.
"""
register!(ctpixels::CtPixels, ct::Int) = 
	pixel -> register!(ctpixels, ct, pixel)

"""
	ispoly(geom) :: Bool

Check if `geom` is a polygon or a multipolygon.
"""
function ispoly(geom)
	GI.isgeometry(geom) || return false
	trait = GI.geomtrait(geom)
	return isa(trait, GI.MultiPolygonTrait) || isa(trait, GI.PolygonTrait)
end

"""
	rasterize(geoms, layer::Raster) :: CtPixels
	rasterize(shptable::AG.IFeatureLayer, layer::Raster) :: CtPixels

For a collection of (multi)polygons, rasterize each of them and write the 
results in a CtPixels.
"""
function rasterize(geoms, layer::Raster)
	all(ispoly.(geoms)) || 
		throw(ArgumentError("`geoms` are not (multi)polygons!"))
	indices = indicate(layer)
	ctpixels = CtPixels(indices)
	for (i, s) = enumerate(geoms)
		rasterize!(indices, s, boundary=:touches, fill=register!(ctpixels,i))
	end
	return ctpixels
end
rasterize(shptable::AG.IFeatureLayer, layer::Raster) = 
	rasterize(AG.getgeom.(shptable), layer)
