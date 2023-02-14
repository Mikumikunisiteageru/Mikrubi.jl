# test/exsim.jl

using Mikrubi
using Test

function county_list(R, hw=240)
	sl = 1
	sr = 1 + R
	pixels = -hw+0.5 : hw-0.5
	junctions = vcat(-hw, reverse(-sl:-sl:1-hw), 0:sr:hw-1+sr) .+ 0.5
	cumsum(in.(pixels, (junctions,))), collect(pixels)
end

@testset "field" begin
	R = 30
	ctids, vars = county_list(R)
	@test length(ctids) == 480
	@test length(vars) == 480
	@test issorted(ctids)
	@test issorted(vars)
	@test isa(ctids, Vector{Int})
	@test sum(ctids) == 87572
	@test isa(vars, Vector{Float64})
	@test sum(vars) == 0.0
	global asym = MikrubiField(ctids, vars, vars)
	@test isa(asym, MikrubiField{Int, Float64, Float64})
	@test asym.mcounty == 248
	@test asym.npixel == 480
	@test asym.dvar == 1
	@test asym.ctids == ctids
	@test asym.locs == repeat(vars, 1, 1)
	@test asym.vars == repeat(vars, 1, 1)
end

@testset "model" begin
	params = [0.02, 0, 1]
	@test_throws DimensionMismatch MikrubiModel(2, params)
	global model = MikrubiModel(1, params)
	@test isa(model, MikrubiModel{Float64})
	@test model.dvar == 1
	@test model.params == params
	pc = probcounties(asym, model)
	@test isapprox(pc[39], 3.253660629809474e-8)
	@test isapprox(pc[239], 0.2687645074317543)
	@test isapprox(pc[248], 1.3446977198405818e-8)
	pp = Mikrubi.probpixels(asym, model)
	@test pp == predict(asym.vars, model)
	@test isapprox(pp[39], 3.253660629809474e-8)
	@test isapprox(pp[239], 0.2687645074317544)
	@test isapprox(pp[248], 0.26454071698645776)
	c = predictcounty(asym, model, 239)
	@test isa(c, Vector{Tuple{Vector{Float64}, Float64}})
	@test length(c) == 1
	@test length(c[1]) == 2
	@test c[1][1] == [-1.5]
	@test isapprox(c[1][2], 0.2687645074317544)
	sample = samplecounties(asym, model)
	@test isa(sample, Vector{Int})
end

@testset "fit" begin
	sample = [183,195,196,203,204,206,207,208,222,233,237,240,241,242]
	optresults = []
	model1 = fit(asym, sample; optresult=optresults)
	optresult = optresults[1]
	@test optresult.ls_success
	@test isapprox(optresult.minimum, 32.021621971313365)
	@test isa(model1, MikrubiModel{Float64})
	@test model1.dvar == 1
	@test all(isapprox.(model1.params, 
		[0.025545912586333833, 0.013468233505695898, 1.1292564512607859]))
end
