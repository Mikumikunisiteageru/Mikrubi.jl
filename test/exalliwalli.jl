# test/exalliwalli.jl

using Mikrubi
using Test

import Rasters

dir = pwd()
cd(joinpath(pkgdir(Mikrubi), "examples", "alliwalli"))

@testset "readfield" begin
	global china = readfield("chinafield.mkuf")
	@test isa(china, MikrubiField{Int, Float64, Float64})
	@test china.mcounty == 2893
	@test china.npixel == 62716
	@test china.dvar == 3
	@test sum(china.ctids) == 84874351
	@test isapprox(sum(china.locs), 8.883390833333332e6)
	@test isapprox(sum(china.vars), -32830.26854044552)
	@test china.ctids[3939] == 199
	@test isapprox(china.locs[3939, :], 
		[105.08333333333331, 25.083333333333343])
	@test isapprox(china.vars[3939, :], 
		[-2.8994803671788847, 0.7973346798555266, -0.5631648131291713])
end

@testset "readlayers" begin
	global ylayers = readlayers("ylayers")
	@test isa(ylayers, Rasters.RasterStack)
	@test size(ylayers) == (2160, 1080, 1)
	@test sum(Rasters.boolmask(ylayers)) == 35813
	ybio1 = ylayers[:ybio1]
	@test isa(ybio1, Rasters.Raster{Float64})
	@test size(ybio1) == (2160, 1080, 1)
	@test sum(Rasters.boolmask(ybio1)) == 35813
	@test isapprox(ybio1[1839, 239, 1], 1.7729962288447991)
	@test isapprox(sum(skipmissing(ybio1)), 1.738049704158584e-10)
end

@testset "readlist" begin
	global ctlist = readlist("countylist.txt")
	@test isa(ctlist, Vector{Int})
	@test length(ctlist) == 46
	@test ctlist[39] == 72
	@test sum(ctlist) == 36844
	ctlistfilename = joinpath(tempdir(), "countylist.txt")
	isfile(ctlistfilename) && rm(ctlistfilename)
	writelist(ctlistfilename, ctlist)
	ctlist1 = readlist(ctlistfilename)
	@test isa(ctlist1, Vector{Int})
	@test ctlist == ctlist1
	rm(ctlistfilename)
end

@testset "fit" begin
	optresults = []
	global model = fit(china, ctlist; optresult=optresults)
	optresult = optresults[1]
	@test optresult.ls_success
	@test isapprox(optresult.minimum, 126.65599400745549)
	@test isa(model, MikrubiModel{Float64})
	@test model.dvar == 3
	@test isapprox(model.params, [1.4842288152354197, 
		-1.3603311815698715, -0.38761691866210646, 1.1231074177981228, 
		 1.2090116395112087, -0.10334796181736790, 14.747024521778938, 
		-14.878922083170924,  11.9705675223002300, 30.299436373642205])
	model1 = MikrubiModel(3, [1.4,-1.4,-0.4,1.1,1.2,-0.1,14.7,-14.9,12.0,30.3])
	e0 = Mikrubi.mlogL(china, ctlist, model.params)
	e1 = Mikrubi.mlogL(china, ctlist, model1.params)
	@test isapprox(e0, 126.65599400745549)
	@test isapprox(e1, 303.59978177848010)
	@test e0 < e1
end

@testset "writemodel" begin
	modelfilename = joinpath(tempdir(), "model.mkum")
	isfile(modelfilename) && rm(modelfilename)
	writemodel(modelfilename, model)
	model1 = readmodel(modelfilename)
	@test model.dvar == model1.dvar
	@test isapprox(model.params, model1.params)
	rm(modelfilename)
end

@testset "predict" begin
	global geodist = predict(ylayers, model)
	@test isa(geodist, Rasters.Raster{Float64})
	@test length(geodist) == 2332800
	@test size(geodist) == (2160, 1080, 1)
	@test sum(Rasters.boolmask(geodist)) == 35813
	@test isapprox(geodist[39*43, 39*10, 1], 0.013895678502362063)
	@test isapprox(sum(skipmissing(geodist)), 28.461626260733837)
end

@testset "writelayer" begin
	layerfilename = joinpath(tempdir(), "geodist.tif")
	isfile(layerfilename) && rm(layerfilename)
	writelayer(layerfilename, geodist)
	geodist1 = Rasters.Raster(layerfilename)
	@test isa(geodist1, Rasters.Raster{Float64})
	@test geodist == geodist1
	rm(layerfilename)
end

@testset "predict" begin
	p = predict(china, model)
	@test isa(p, Dict{Int, Vector{Tuple{Vector{Float64}, Float64}}})
	@test length(p) == 2893
	@test sum(length.(values(p))) == 62716
	p39 = p[39]
	@test isa(p39, Vector{Tuple{Vector{Float64}, Float64}})
	@test length(p39) == 30
	loc, prob = first(p39)
	@test isapprox(loc, [100.08333333333331, 28.583333333333336])
	@test isapprox(prob, 0.014382846487993373)
end

@testset "probcounties" begin
	pc = probcounties(china, model)
	@test isa(pc, Dict{Int, Float64})
	@test length(pc) == 2893
	@test sum(keys(pc)) == 4188157
	@test isapprox(sum(values(pc)), 45.27843660370468)
	@test isapprox(pc[39], 0.19286274249159574) 
end

@testset "predictcounty" begin
	p39 = predictcounty(china, model, 39)
	@test isa(p39, Vector{Tuple{Vector{Float64}, Float64}})
	@test length(p39) == 30
	loc, prob = first(p39)
	@test isapprox(loc, [100.08333333333331, 28.583333333333336])
	@test isapprox(prob, 0.014382846487993373)
	@test p39 == predict(china, model)[39]
end

cd(dir)
