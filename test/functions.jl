# test/functions.jl

using Mikrubi
using Test

@testset "textwrap and @tw_str" begin
	@test Mikrubi.textwrap("Very   \n\t   good") == "Very good"
	@test Mikrubi.textwrap("Very 
		good") == "Very good"
	@test Mikrubi.textwrap("Very
		good") == "Very good"
	@test Mikrubi.tw"Very 
		good" == "Very good"
	@test Mikrubi.tw"Very
		good" == "Very good"
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
	colmean, projwstd = Mikrubi.princompvars(matrix)
	pcamatrix = (matrix .- colmean) * projwstd
	@test all(isapprox.((pcamatrix' * pcamatrix)[[2,3,4,6,7,8]], 0, atol=1e-10))
	@test isapprox(sum(pcamatrix), 0, atol=1e-10)
	@test isapprox(sum(pcamatrix.^2), 12, atol=1e-10)
end

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
