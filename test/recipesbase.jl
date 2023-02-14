# test/recipesbase.jl

using Mikrubi
using Test

using Plots
# using FileIO

import GADM
import RasterDataSources; const RDS = RasterDataSources

shppath = GADM.download("NPL")

get!(ENV, "RASTERDATASOURCES_PATH", tempdir())
RDS.getraster(RDS.WorldClim{RDS.BioClim}, res="10m")
climpath = RDS.rasterpath(RDS.WorldClim{RDS.BioClim})

shptable = readshape(shppath, 1)
layers = readlayers(climpath)
layer = first(layers)
ctpixels = rasterize(shptable, layer)
field, ylayers = makefield(layers, shptable)

white = RGB(1.0, 1.0, 1.0)

# function checksumifrepl(filename, r, g, b)
	# savefig(joinpath(tempdir(), filename))
	# img = FileIO.load(joinpath(tempdir(), filename))
	# sumcolor = sum(white .- img)
	# @test isapprox(sumcolor, RGB(r, g, b))
# end

@testset "plot a shptable" begin
	@test isa(plot(shptable), Plots.Plot)
	# checksumifrepl("shptable.png", 
		# 41932.725490196084, 39602.3019607843, 41112.08627450979)
	@test isa(plot!(shptable), Plots.Plot)
end

@testset "plot a layer" begin
	@test isa(plot(layer), Plots.Plot)
	# checksumifrepl("layer.png", 
		# 18255.329411764702, 22052.23921568627, 26050.831372549022)
	@test isa(plot!(layer), Plots.Plot)
end

@testset "plot a CtPixels" begin
	@test isa(plot(ctpixels), Plots.Plot)
	# checksumifrepl("ctpixels.png", 
		# 26157.454901960777, 24958.325490196075, 23180.349019607835)
	@test isa(plot!(ctpixels), Plots.Plot)
	@test isa(plot(ctpixels; salt=21), Plots.Plot)
	# checksumifrepl("ctpixels2.png", 
		# 25980.529411764706, 25731.533333333326, 22352.149019607838)
	@test isa(plot!(ctpixels; salt=21), Plots.Plot)
end

@testset "plot a Mikrubi field" begin
	@test isa(plot(layer, field), Plots.Plot)
	# checksumifrepl("field.png", 
		# 29293.13725490196, 27610.82745098039, 30779.545098039212)
	@test isa(plot!(layer, field), Plots.Plot)
	@test isa(plot(layer, field; func=identity), Plots.Plot)
	# checksumifrepl("field2.png", 
		# 25963.666666666664, 30761.10588235294, 31637.52549019608)
	@test isa(plot!(layer, field; func=identity), Plots.Plot)
end

Plots.closeall()
