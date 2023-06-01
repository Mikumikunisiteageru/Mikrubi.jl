# test/exjui.jl

using Mikrubi
using Logistics
using Test

import ArchGDAL; const AG = ArchGDAL
import GADM
import GDAL
import RasterDataSources; const RDS = RasterDataSources
import Rasters

@testset "GADM" begin
	global shppath = GADM.download("NPL"; version="4.1")
	@test isdir(shppath)
	shpfiles = readdir(shppath, join=false, sort=true)
	@test length(shpfiles) == 1
	global gpkg, = shpfiles
	@show shppath
	@show gpkg
	@test gpkg == "gadm41_NPL.gpkg"
end

@testset "RasterDataSources" begin
	get!(ENV, "RASTERDATASOURCES_PATH", tempdir())
	RDS.getraster(RDS.WorldClim{RDS.BioClim}, res="10m")
	global climpath = RDS.rasterpath(RDS.WorldClim{RDS.BioClim})
	@test isdir(climpath)
	filenames = readdir(climpath, join=true)
	@test length(filenames) == 20
	@test length(filter(isfile, filenames)) == 19
	@test length(filter(isdir, filenames)) == 1
	@test unique(last.(splitext.(filter(isfile, filenames)))) == [".tif"]
	@test last(splitpath(first(filenames))) == "wc2.1_10m_bio_1.tif"
end

@testset "readshape" begin
	@test_throws ErrorException readshape(joinpath(shppath, gpkg))
	@test_throws ErrorException readshape(joinpath(shppath, gpkg), -1)
	@test isa(readshape(joinpath(shppath, gpkg), 4), AG.IFeatureLayer)
	@test_throws ErrorException readshape(joinpath(shppath, gpkg), 5)
	@test_throws ErrorException readshape(shppath; extset=[])
	@test_throws ErrorException readshape(shppath)
	@test_throws ErrorException readshape(shppath, -1)
	body, tail = splitext(gpkg)
	gpkg_ = body * "_" * tail
	cp(joinpath(shppath, gpkg), joinpath(shppath, gpkg_))
	msg = Mikrubi.tw"Multiple files with extensions in `extset` exist 
		in the directory! Now choose a random one to read."
	@test_logs (:warn, msg) readshape(shppath, 3)
	rm(joinpath(shppath, gpkg_))
	global shptable = readshape(shppath, 3)
	@test isa(shptable, AG.IFeatureLayer)
	@test AG.getname(shptable) == "ADM_ADM_3"
	@test AG.nfeature(shptable) == 75
	pt = AG.getgeom(AG.getgeom(AG.getgeom(first(shptable), 0), 0), 0)
	@test isapprox(AG.getx(pt, 0), 85.406402588)
	@test isapprox(AG.gety(pt, 0), 27.632347107)
end

@testset "goodcolumns" begin
	gc = goodcolumns(shptable)
	@test sort!(collect(keys(gc))) == ["GID_3", "NAME_3"]
end

@testset "lookup" begin
	@test_throws ErrorException lookup(shptable)
	@test_throws ErrorException lookup(shptable, "Any", 1)
	@test_throws ErrorException lookup(shptable, :Any, 1)
	@test nothing === lookup(shptable, "GID_3", 1)
end

@testset "readlayers" begin
	global layers = readlayers(climpath)
	@test isa(layers, Rasters.RasterStack)
	@test length(layers) == 19
	@test keys(layers)[1] == Symbol("wc2.1_10m_bio_1")
	@test keys(layers)[19] == Symbol("wc2.1_10m_bio_19")
	bm = Rasters.boolmask(layers)
	@test sum(bm) == 808053
	@test length(collect(skipmissing(first(layers)))) == 808053
	@test isapprox(sum(collect(skipmissing(first(layers)))), -3.262779f6)
	@test isapprox(collect(skipmissing(first(layers)))[1:10], Float32[0.0, 0.0, 
		-2.5923913, -8.346475, -16.416666, -17.895636, 
		-6.286264, -6.5873747, -5.2716227, -2.6854463])
end

@testset "rasterize" begin
	layer = first(layers)
	@test_throws MethodError Mikrubi.CtPixels()
	@test_throws MethodError Mikrubi.CtPixels(layer)
	indices = Mikrubi.indicate(layer)
	@test isa(indices, Rasters.Raster{Int})
	@test indices[0393939] == 0
	@test indices[1393939] == 1393939
	cpx = Mikrubi.CtPixels(indices)
	@test isa(cpx, Mikrubi.CtPixels)
	@test length(cpx) == 0
	global ctpixels = rasterize(shptable, layer)
	@test isa(ctpixels, Mikrubi.CtPixels)
	@test length(ctpixels) == 1127
	@test ctpixels.indices == Mikrubi.indicate(layer) == cpx.indices
	@test length(ctpixels.list) == 1127
	@test ctpixels.list[39] == (5, 811592)
	@test Mikrubi.getcounty(ctpixels, 39) == 5
	@test Mikrubi.getpixel(ctpixels, 39) == 811592
	@test Mikrubi.getcounties(ctpixels)[39] == 5
	@test Mikrubi.getpixels(ctpixels)[39] == 811592
	@test isa(Mikrubi.getcounties(ctpixels), Vector{Int})
	@test isa(Mikrubi.getpixels(ctpixels), Vector{Int})
	@test length(Mikrubi.getcounties(ctpixels)) == 1127
	@test length(Mikrubi.getpixels(ctpixels)) == 1127
end

@testset "buildfield" begin
	layers_ = deepcopy(layers)
	ctpixels_ = deepcopy(ctpixels)
	ctpixels_.indices[1, 1, 1] = 1
	push!(ctpixels_.list, (76, 1))
	Mikrubi.masklayers!(layers_, ctpixels_)
	dimlower = Mikrubi.DimLower(layers_)
	idx, ematrix, elayers = dimlower(layers_; rabsthres=0.8, nprincomp=3)
	field = Mikrubi.buildfield(ctpixels_, idx, ematrix, first(layers_))
	@test isa(field, MikrubiField)
end

@testset "DimLower" begin
	layers_ = deepcopy(layers)
	Mikrubi.masklayers!(layers_, ctpixels)
	matrix, idx = Mikrubi.extractlayers(layers_)
	@test matrix[:, 1] == collect(skipmissing(first(layers_)))
	colid = Mikrubi.selectvars(matrix; rabsthres=0.8)
	@test colid == [2, 3, 5, 7, 12, 14, 17, 19]
	smatrix = matrix[:, colid]
	@test isapprox(smatrix[1:10], Float32[11.106563, 12.474979, 
		11.979521, 12.020729, 12.416375, 10.518167, 
		11.936396, 11.022813, 11.316313, 11.497708])
	colmean, projwstd = Mikrubi.princompvars(smatrix; nprincomp=3)
	@test isapprox(colmean[:], Float32[11.103582, 44.19512, 
		24.790443, 25.243826, 1304.4342, 7.5794873, 47.97265, 81.61539])
	@test isapprox(projwstd, Float32[
		0.4013147 0.3629837 0.18809055; 
		-0.061730698 -0.13327469 -0.108475216; 
		0.030078402 -0.036740705 0.058986623; 
		0.13204093 0.17373833 0.11315063; 
		0.00021478014 -0.00088945474 0.00058913097; 
		-0.13801464 -0.013789057 0.047764596; 
		-0.01994192 0.0015393574 0.021920582; 
		-0.008666443 0.008949931 0.0048636016])
	f_, yl_, dimlower = Mikrubi._makefield(layers, ctpixels, 0.8, 3)
	@test isa(dimlower, Mikrubi.DimLower{Float64})
	@test dimlower.new == false
	@test dimlower.colid == [2, 3, 5, 7, 12, 14, 17, 19]
	@test isa(dimlower.colmean, Matrix{Float64})
	@test size(dimlower.colmean) == (1, 8)
	@test isapprox(dimlower.colmean[:], Float32[11.103582, 44.19512, 
		24.790443, 25.243826, 1304.4342, 7.5794873, 47.97265, 81.61539])
	@test isa(dimlower.projwstd, Matrix{Float64})
	@test size(dimlower.projwstd) == (8, 3)
	@test isapprox(dimlower.projwstd, Float32[
		0.4013147 0.3629837 0.18809055; 
		-0.061730698 -0.13327469 -0.108475216; 
		0.030078402 -0.036740705 0.058986623; 
		0.13204093 0.17373833 0.11315063; 
		0.00021478014 -0.00088945474 0.00058913097; 
		-0.13801464 -0.013789057 0.047764596; 
		-0.01994192 0.0015393574 0.021920582; 
		-0.008666443 0.008949931 0.0048636016])
	idx, ematrix, elayers = dimlower(layers; rabsthres=0.8, nprincomp=3)
	@test isa(idx, Vector{Int})
	@test length(idx) == 808053
	@test sum(idx) == 923357905365
	@test isa(ematrix, Matrix{Float64})
	@test size(ematrix) == (808053, 3)
	@test isapprox(sum(ematrix), 3.4257468f6)
end

@testset "makefield" begin
	f0, yl0 = makefield(layers, ctpixels)
	@test sprint(show, f0) === 
		"Mikrubi Field: geo_dim = 2, env_dim = 3, 1127 pixels, and 75 counties"
	@test sum(Rasters.boolmask(layers)) == 808053
	@test isa(f0, MikrubiField)
	@test isa(yl0, Rasters.RasterStack)
	path = joinpath(tempdir(), "field.mkuf")
	@test nothing === writefield(path, f0)
	str = String(read(path))
	@test str[1:12] == "I\tL\tL\tV\tV\tV\n"
	pathformula(s) = joinpath(tempdir(), "layer_$s.tif")
	@test nothing === writelayers(pathformula.(5:7), yl0)
	yl0_ = readlayers(pathformula.(5:7))
	@test length(yl0_) == 3
	@test size(yl0_) == (2160, 1080, 1)
	@test nothing === writelayers(pathformula("*"), yl0)
	yl0_ = readlayers(pathformula.(1:3))
	@test length(yl0_) == 3
	@test size(yl0_) == (2160, 1080, 1)
	global field, ylayers = makefield(layers, shptable);
	@test sum(Rasters.boolmask(layers)) == 808053
	@test yl0 == ylayers
	@test f0.ctids == field.ctids
	@test f0.locs == field.locs
	@test f0.vars == field.vars
	@test sum(field.ctids) == 45008
	@test isapprox(sum(field.locs), 126492.66666666664)
	@test isapprox(sum(field.vars), -395.2514f0)
	@test length(ylayers) == 3
	bm = Rasters.boolmask(ylayers)
	@test sum(bm) == 585
	@test isapprox(collect(maximum(ylayers)), 
		Float32[3.5854542, 4.006934, 2.3392994])
	@test isapprox(sum.(collect(ylayers[bm])), zeros(3), atol=1e-10)
	f1, yl1, pyl1 = makefield(layers, ctpixels, layers)
	@test sum(Rasters.boolmask(layers)) == 808053
	@test yl0 == yl1
	@test isapprox(collect(yl1[bm]), collect(pyl1[bm]))
	@test sum(Rasters.boolmask(pyl1)) == 808053
	f2, yl2, pyl2 = makefield(layers, shptable, layers)
	@test sum(Rasters.boolmask(layers)) == 808053
	@test yl0 == yl2
	@test isapprox(collect(yl2[bm]), collect(pyl2[bm]))
	@test sum(Rasters.boolmask(pyl2)) == 808053
end

@testset "lookup" begin
	regions = ["Dadeldhura", "Doti", "Bajhang", "Kalikot", 
		"Mugu", "Jajarkot", "Jumla", "Rolpa", "Dolpa", "Baglung", 
		"Mustang", "Manang", "Gorkha", "Nuwakot", "Rasuwa", 
		"Okhaldhunga", "Solukhumbu"]
	@test length(regions) == 17
	@test allunique(regions)
	@test lookup(shptable, "NAME_3", "Bajhang") == 41
	@test AG.getfeature(x -> AG.getfield(x, 9), shptable, 41) == "Bajhang"
	@test isnothing(lookup(shptable, "NAME_3", "Bazhang"))
	@test lookup(shptable, "NAME_2", "Seti") == [40, 41, 42, 43, 44]
	global regcodes = lookup(shptable, "NAME_3", regions)
	@test isa(regcodes, Vector{Int})
	@test regcodes == [37,43,41,53,54,48,52,57,50,60,61,67,64,6,7,31,34]
end

@testset "findnearest" begin
	@test_throws ErrorException Mikrubi.findnearest([85.32], field)
	@test_throws ErrorException Mikrubi.findnearest([8, 5, 3, 2], field)
	@test Mikrubi.findnearest([85.321201, 27.722903], field) == 19
	@test Mikrubi.findnearest([85.321201 27.722903], field) == 19
	@test Mikrubi.findnearest([85.321201 27.722903]', field) == 19
	@test isapprox(field.locs[19, :], [85.25, 27.75])
end

@testset "findnearests" begin
	@test Mikrubi.findnearests(
		[85.321201 27.722903; 89.368606 24.85836], field) == [19, 345]
	@test Mikrubi.findnearests(
		[[85.321201, 27.722903], [89.368606, 24.85836]], field) == [19, 345]
end

@testset "fit" begin
	optresults = []
	global model = fit(field, regcodes; optresult=optresults)
	optresult = optresults[1]
	@test optresult.ls_success
	@test isapprox(optresult.minimum, 20.805983368146116)
	@test isa(model, MikrubiModel{Float64})
	@test model.dvar == 3
	@test isapprox(model.params, [7.180739724704129, 
		-0.04075956831789931, -0.54476207315363, 0.8879548516254412, 
		0.3960510254962835, 0.00011517895691697269, 193.3943652333833, 
		-1.5361309085483503, 24.225666096714022, 181.11673123077227])
	model1 = MikrubiModel(3, model.params .* 1.01f0)
	e0 = Mikrubi.mlogL(field, regcodes, model.params)
	e1 = Mikrubi.mlogL(field, regcodes, model1.params)
	@test isapprox(e0, 20.805983368146116)
	@test isapprox(e1, 34.81738959663534)
	@test e0 < e1
	pp = Mikrubi.ppresence(field, model.params)
	@test isa(pp, Vector{Logistic{Float64}})
	@test size(pp) == (1127,)
	pp1 = Mikrubi.ppresence(field.vars, model.params)
	@test pp1 == pp
	@test typeof(pp1) == typeof(pp)
end

@testset "predict" begin
	global geodist = predict(ylayers, model)
	@test isa(geodist, Rasters.Raster{Float64})
	@test length(geodist) == 2332800
	@test size(geodist) == (2160, 1080, 1)
	@test sum(Rasters.boolmask(geodist)) == 585
	@test isapprox(geodist[809439], 0.0016682358794606333)
	@test isapprox(sum(skipmissing(geodist)), 12.929364196004789)
end
