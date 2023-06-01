# test/pyplot.jl

using Mikrubi
using Test

@testset "pyplot functions undefined" begin
	@test ! @isdefined showlayer
	@test ! @isdefined showfield
	@test ! @isdefined showctpixels
	@test ! @isdefined showshptable
end

using PyPlot
using PyCall

@testset "pyplot functions defined" begin
	@test @isdefined showlayer
	@test @isdefined showfield
	@test @isdefined showctpixels
	@test @isdefined showshptable
end

import GADM
import RasterDataSources; const RDS = RasterDataSources

import ArchGDAL; const AG = ArchGDAL

shppath = GADM.download("NPL")

get!(ENV, "RASTERDATASOURCES_PATH", tempdir())
RDS.getraster(RDS.WorldClim{RDS.BioClim}, res="10m")
climpath = RDS.rasterpath(RDS.WorldClim{RDS.BioClim})

shptable = readshape(shppath, 3)
layers = readlayers(climpath)
layer = first(layers)
ctpixels = rasterize(shptable, layer)
field, ylayers = makefield(layers, shptable)

@testset "setplot" begin
	@test setplot(PyPlot) === nothing
end

@testset "geom2mat" begin
	multipolygon = AG.getgeom(first(shptable), 0)
	@test isa(multipolygon, AG.IGeometry{AG.wkbMultiPolygon})
	parts, mat = Mikrubi.geom2mat(multipolygon)
	@test parts == [0]
	@test size(mat) == (2, 118)
	@test isapprox(last(mat), 27.632347107)
	polygon = AG.getgeom(multipolygon, 0)
	@test isa(polygon, AG.IGeometry{AG.wkbPolygon})
	parts, mat = Mikrubi.geom2mat(polygon)
	@test parts == [0]
	@test size(mat) == (2, 118)
	@test isapprox(last(mat), 27.632347107)
	linearring = AG.getgeom(polygon, 0)
	@test isa(linearring, AG.IGeometry{AG.wkbLineString})
	@test_throws ErrorException Mikrubi.geom2mat(linearring)
end

@testset "showshptable" begin
	@test isa(showshptable(shptable), Vector{<:PyObject})
	close()
end

@testset "showlayer" begin
	@test isa(showlayer(layer), PyObject)
	close()
end

@testset "showctpixels" begin
	@test isa(showctpixels(ctpixels), PyObject)
	close()
	@test isa(showctpixels(ctpixels, layer), PyObject)
	close()
end

@testset "showfield" begin
	@test isa(showfield(field), PyObject)
	close()
	@test isa(showfield(field, layer), PyObject)
	close()
end
