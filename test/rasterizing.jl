# src/rasterizing.jl

using Mikrubi
using Test

import ArchGDAL; const AG = ArchGDAL
using Rasters
using .LookupArrays

sintv(n) = Sampled(1:n; sampling=Intervals(End()))

@testset "allium-shape" begin
	raster = Raster(zeros(X(sintv(7)), Y(sintv(9))))
	@test isapprox([Mikrubi.centercoords(raster, 39)...], [3.5, 5.5])
	indices = Mikrubi.indicate(raster)
	wkt = "POLYGON ((0.5 0.5, 2.2 5.3, 1 8, 2.3 8.7, " * 
		"3.5 6.5, 5.5 8.5, 6.5 7.5, 4.2 4.8, 2.8 0, 0.5 0.5))"
	geom = AG.fromWKT(wkt)
	ctpixels = Mikrubi.rasterize([geom], raster)
	@test isa(ctpixels, Mikrubi.CtPixels)
	@test all(isone.(Mikrubi.getcounties(ctpixels)))
	pixels = Mikrubi.getpixels(ctpixels)
	body = vcat(1:4, 8:11, 16:18, 23:25, 30:33, 37:41, 44:48, 51:56, 58:59, 62)
	skin = vcat(56, 49, 61, 63)
	@test isempty(setdiff(body, pixels))
	@test isempty(setdiff(pixels, union(body, skin)))
end

@testset "ribbon-shape" begin
	raster = Raster(zeros(X(sintv(8)), Y(sintv(8))))
	@test isapprox([Mikrubi.centercoords(raster, 39)...], [6.5, 4.5])
	indices = Mikrubi.indicate(raster)
	# wkt = Mikrubi.tw"MULTIPOLYGON
		# (((1 4, 3.5 7, 3.9 7, 1.4 4, 3.9 1, 3.5 1, 1 4)), 
		 # ((1.6 4, 4.1 7, 4.5 7, 7 4, 4.5 1, 4.1 1, 1.6 4), 
		  # (2 4, 4 1.6, 6 4, 4 6.4, 2 4), 
		  # (6.4 4, 4.2 1.36, 4.3 1.24, 6.6 4, 4.3 6.76, 4.2 6.64, 6.4 4)))"
	wkt = Mikrubi.tw"MULTIPOLYGON
		(((1.01 4.01, 3.5 7, 3.9 7, 1.4 4, 3.9 1, 3.5 1.01, 1.01 4.01)), 
		 ((1.6 4.01, 4.1 7, 4.5 7, 7 4, 4.5 1, 4.1 1.01, 1.6 4.01), 
		  (2.01 4.01, 4.01 1.6, 6 4, 4 6.4, 2.01 4.01), 
		  (6.4 4, 4.2 1.36, 4.3 1.24, 6.6 4, 4.3 6.76, 4.2 6.64, 6.4 4)))"
	# A bug in `_fill_line!` of Rasters v0.5.1, see: 
	# https://github.com/rafaqz/Rasters.jl/issues/376
	geoms = AG.fromWKT(wkt)
	ctpixels = Mikrubi.rasterize([geoms], raster)
	@test isa(ctpixels, Mikrubi.CtPixels)
	@test all(isone.(Mikrubi.getcounties(ctpixels)))
	pixels = Mikrubi.getpixels(ctpixels)
	body = vcat(11:14, 18:23, 26:27, 30:31, 34:35, 38:39, 42:47, 51:54)
	skin = vcat(4, 5, 25, 32, 33, 40, 60, 61)
	@test isempty(setdiff(body, pixels))
	@test isempty(setdiff(pixels, union(body, skin)))
	geom1 = AG.getgeom(geoms, 0)
	geom2 = AG.getgeom(geoms, 1)
	ctpixels = Mikrubi.rasterize([geom1, geom2], raster)
	@test isa(ctpixels, Mikrubi.CtPixels)
	@test sort!(unique(Mikrubi.getcounties(ctpixels))) == [1, 2]
	pixels1 = last.(filter(x -> x[1] == 1, ctpixels.list))
	body1 = vcat(11:12, 18:20, 26:27, 34:35, 42:44, 51:52)
	skin1 = vcat(4, 25, 33, 60)
	@test isempty(setdiff(body1, pixels1))
	@test isempty(setdiff(pixels1, union(body1, skin1)))
	pixels2 = last.(filter(x -> x[1] == 2, ctpixels.list))
	body2 = vcat(12:14, 19:23, 26:27, 30:31, 34:35, 38:39, 43:47, 52:54)
	skin2 = vcat(5, 32, 40, 61)
	@test isempty(setdiff(body2, pixels2))
	@test isempty(setdiff(pixels2, union(body2, skin2)))
end
