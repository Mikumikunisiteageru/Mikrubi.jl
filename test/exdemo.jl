# test/exdemo.jl

using Mikrubi
using Test

@testset "a small example" begin
	cty = ["a", "b", "c", "d", "d", "e", "f", "f", "e"]
	geo = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
	env = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
	field = MikrubiField(cty, geo, env)
	@test field.npixel  == 9
	@test field.mcounty == 6
	@test field.dvar    == 2
	@test field.ids == string.(collect("abcdef"))
	@test getindex.([field.starts], field.ids) == [1, 2, 3, 4, 6, 8]
	@test getindex.([field.stops],  field.ids) == [1, 2, 3, 5, 7, 9]
	model = fit(field, ["a", "c", "f"])
	@test model.dvar == 2
	@test 0 <= length(samplecounties(field, model)) <= 6
	@test isapprox(sum(Mikrubi.probpixels(field, model)), 3.137251f0)
	@test isapprox(sum(values(probcounties(field, model))), 2.8943691f0)
	@test isapprox(lipschitz(model, field, wholespace=false), 0.21698886f0)
	@test isapprox(lipschitz(model, field, wholespace=true), 4.4162006f0)
end
