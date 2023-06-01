# test/pyplot.jl

using Mikrubi
using Test

using PyPlot
using PyCall
# test this manually because I don't want to include them in Project.toml

import GADM
import RasterDataSources; const RDS = RasterDataSources

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
	@test_throws ErrorException showshptable(shptable)
	@test setplot(PyPlot) == PyPlot
end

@testset "show a shptable" begin
	@test isa(showshptable(shptable), Vector{<:PyObject})
	close()
end

@testset "show a layer" begin
	@test isa(showlayer(layer), PyObject)
	close()
end

@testset "show a CtPixels" begin
	@test isa(showctpixels(ctpixels), PyObject)
	@test isa(showctpixels(ctpixels, layer), PyObject)
	close()
end

@testset "show a Mikrubi field" begin
	@test isa(showfield(field), PyObject)
	@test isa(showfield(field, layer), PyObject)
	close()
end
