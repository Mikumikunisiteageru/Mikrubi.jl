# test/runtests.jl

using Mikrubi
using Test
import ArchGDAL
const AG = ArchGDAL
using GeoArrays

J(args...) = Mikrubi.I.(args...)

@testset "PROGRAMMING TRICKS AND GADGETS" begin

	@testset "I" begin
		@test Mikrubi.I.(1, [2,3,4]) == [(1,2), (1,3), (1,4)]
		@test Mikrubi.I.([1,0,-1], [2,3,4]) == [(1,2), (0,3), (-1,4)]
	end

	@testset "textwrap" begin
		@test Mikrubi.textwrap("Very   \n\t   good") == "Very good"
	end

	@testset "allsame" begin
		@test Mikrubi.allsame([1, 1, 2]) == false
		@test Mikrubi.allsame([1, 1, 1]) == true
		@test Mikrubi.allsame([1]) == true
		@test_throws BoundsError Mikrubi.allsame([])
	end

	@testset "colmatrix" begin
		@test Mikrubi.colmatrix([1]) == fill(1, 1, 1)
		@test Mikrubi.colmatrix([1, 1]) == fill(1, 2, 1)
		@test Mikrubi.colmatrix([1 1; 1 1]) == fill(1, 2, 2)
	end

end

@testset "MATHEMATIC TOOLS" begin

	@testset "logistic" begin
		@test logistic(-Inf) == 0.0
		@test logistic(+Inf) == 1.0
		@test logistic(0.0) == 0.5
		@test logistic(0.39) + logistic(-0.39) == 1.0
	end

	@testset "loglogistic" begin
		@test loglogistic(-Inf) == -Inf
		@test loglogistic(+Inf) == 0.0
		@test loglogistic(0.0) == -log(2.0)
	end

end

@testset "FUNCTIONS FOR RASTERIZING" begin

	@testset "interpolate" begin
		interpolate = Mikrubi.interpolate
		@test interpolate(1.0, 2.0, 3.0, 4.0)(0.0) == 2.0
		@test interpolate(1, 2, 3, 4)(0) == 2.0
		@test interpolate(1, 2, 3, 3)(0) == 3.0
		@test isinf(interpolate(1, 1, 3, 4)(0))
		@test isinf(interpolate(1, 1, 3, 4)(2))
	end

	@testset "RASTERIZING A POLYGON OF Y-SHAPE" begin
		wkt = "POLYGON ((0.5 0.5, 2.2 5.3, 1 8, 2.3 8.7, " *
			"3.5 6.5, 5.5 8.5, 6.5 7.5, 4.2 4.8, 2.8 0, 0.5 0.5))"
		geom = AG.fromWKT(wkt)
		parts, mat = Mikrubi.geom2mat(geom)
		layer = GeoArray(zeros(Union{Missing, Float32}, 7, 9), 
			GeoArrays.AffineMap(Float32[1 0; 0 1], Float32[0, 0]))
			# GeoArrays.AffineMap(Float32[0 1; 1 0], Float32[0, 0]))
		@test parts == [0]
		@test mat == [0.5 2.2 1.0 2.3 3.5 5.5 6.5 4.2 2.8 0.5; 
					  0.5 5.3 8.0 8.7 6.5 8.5 7.5 4.8 0.0 0.5]
		@test sort(Mikrubi.scanline(mat, parts.+1, 7, 9), by=last) ==
			vcat(J.([1:3, 2:3, 2:4, 3:4, 3:4, 3:5, [3,5,6], [2,3,5,6]], 1:8)...)	
		@test sort(collect(Mikrubi.strokepath(mat, parts.+1, 7, 9))) ==
			vcat(J.(1:7, 
			[1:2, 1:9, [1,5,6,8,9], [1:5...,7,8], [5,6,8], [6,7,9], [7,8]])...)
		@test sort(rasterize(geom, layer)) == 
			vcat(J.(1:7, [1:2, 1:9, 1:9, 1:8, 5:8, 6:9, 7:8])...)
	end

	@testset "RASTERIZING A MULTIPOLYGON WITH TWO HOLES" begin
		wkt = "MULTIPOLYGON
			(((1 4, 3.5 7, 3.9 7, 1.4 4, 3.9 1, 3.5 1, 1 4)), 
			 ((1.6 4, 4.1 7, 4.5 7, 7 4, 4.5 1, 4.1 1, 1.6 4), 
			  (2 4, 4 1.6, 6 4, 4 6.4, 2 4), 
			  (6.4 4, 4.2 1.36, 4.3 1.24, 6.6 4, 4.3 6.76, 4.2 6.64, 6.4 4)))"
		geom = AG.fromWKT(wkt)
		parts, mat = Mikrubi.geom2mat(geom)
		layer = GeoArray(zeros(Union{Missing, Float32}, 8, 8), 
			GeoArrays.AffineMap(Float32[1 0; 0 1], Float32[0, 0]))
		@test parts == [0, 7, 14, 19]
		@test mat == 
			[1.0 3.5 3.9 1.4 3.9 3.5 1.0 1.6 4.1 4.5 7.0 4.5 4.1 0+
			 1.6 2.0 4.0 6.0 4.0 2.0 6.4 4.2 4.3 6.6 4.3 4.2 6.4; 
			 4.0 7.0 7.0 4.0 1.0 1.0 4.0 4.0 7.0 7.0 4.0 1.0 1.0 0+
			 4.0 4.0 1.6 4.0 6.4 4.0 4.0 1.36 1.24 4.0 6.76 6.64 4.0]
		@test sort(Mikrubi.scanline(mat, parts.+1, 8, 8)) ==
			[(2,4), (2,5), (3,3), (3,6), (6,3), (6,6), (7,4), (7,5)]
		@test sort(collect(Mikrubi.strokepath(mat, parts.+1, 8, 8))) ==
			vcat(J.(2:7, [3:6, 2:7, [2,3,6,7], [2,3,6,7], 2:7, 3:6])...)
		@test sort(rasterize(geom, layer)) == 
			vcat(J.(2:7, [3:6, 2:7, [2,3,6,7], [2,3,6,7], 2:7, 3:6])...)
	end

	@testset "ABOUT CtPixels" begin
		ctpixels = [((0,1),2), ((3,4),5)]
		@test isa(ctpixels, Mikrubi.CtPixels)
		@test Mikrubi.getpixels(ctpixels) == [(0,1), (3,4)]
		@test Mikrubi.getcounties(ctpixels) == [2, 5]
	end

end

@testset "FUNCTIONS FOR HANDLING RASTER LAYERS" begin

	@testset "dftraverse!" begin
		beststate = Vector(undef, 1);
		bestscore = [(0, 0.0)];
		Mikrubi.dftraverse!(beststate, bestscore, Int[], (0, 0.0), 1, 3,
			Bool[0 0 1; 0 0 0; 1 0 0], [0.0 0.6 0.3; 0.6 0.0 0.9; 0.3 0.9 0.0])
		@test beststate == [[1, 2]]
		@test bestscore == [(2, -0.6)]
	end
	
	@testset "selectvars" begin
		@test Mikrubi.selectvars(
			[1. 4. 7.; 2. 5. 8.; 3. 9. 27.], rabsthres=0.9) == [1, 3]
	end
	
	@testset "princompvars" begin
		matrix = Float64[1 2 3 4; 0 1 2 3; 0 0 1 2; 0 0 0 1]
		pcamatrix, (colmean, projwstd) = Mikrubi.princompvars(matrix)
		@test isapprox(pcamatrix, (matrix .- colmean) * projwstd)
		@test all(isapprox.((pcamatrix' * pcamatrix)[[2,3,4,6,7,8]], 0, atol=1e-10))
		@test isapprox(sum(pcamatrix), 0, atol=1e-10)
		@test isapprox(sum(pcamatrix.^2), 12, atol=1e-10)
	end
	
end

@testset "THE CORE OF MIKRUBI" begin

	@testset "dvar2dparam" begin
		@test Mikrubi.dvar2dparam(1) == 3
		@test Mikrubi.dvar2dparam(2) == 6
		@test Mikrubi.dvar2dparam(3) == 10
		@test Mikrubi.dvar2dparam(19) == 210
	end

	@testset "decomparams" begin
		@test Mikrubi.decomparams([1,2,3], 1) == (fill(1, 1, 1), [2], 3)
		@test Mikrubi.decomparams([1,2,3,4,5,6], 2) == ([1 0; 2 3], [4, 5], 6)
		@test Mikrubi.decomparams([1,2,3,4,5,6,7,8,9,10], 3) == 
			([1 0 0; 2 3 0; 4 5 6], [7, 8, 9], 10)
	end
	
	@testset "A SMALL EXAMPLE" begin
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

end
