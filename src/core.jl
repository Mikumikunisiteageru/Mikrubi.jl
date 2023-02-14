# src/core.jl

"""
	logistic(x)

Compute ``\\operatorname{logistic}(x) := 1 / (1 + \\exp(x))``.
"""
logistic(x::T) where {T<:AbstractFloat} = one(T) / (one(T) + exp(-x))
logistic(x::Real) = logistic(float(x))

"""
	loglogistic(x)

Compute ``\\log(\\operatorname{logistic}(x)) = -\\log(1 + \\exp(-x))``.
"""
loglogistic(x::Real) = -log1p(exp(-x))

"""
	dvar2dparam(dvar::Int) :: Int

Convert dimensionality of an environmental space to the dimensionality of
the induced parameter space, i.e., compute the degrees of freedom for 
positive-definite quadratic functions mapping a `dvar`-dimensional linear space
into real numbers.

# Examples
```julia
julia> dvar2dparam(1)
3

julia> dvar2dparam(3)
10
```
"""
dvar2dparam(dvar::Int) = ((dvar + 1) * (dvar + 2)) >> 1

"""
	MikrubiField{T, U <: Real, V <: AbstractFloat}

	MikrubiField(ctids, locs, vars)

Construct a Mikrubi field containing a number of pixels or points, using the
arguments which should always have the same number of rows
- `ctids::Vector`: a vector containing the county identifiers
- `locs::Array{<:Real}`: an array of geographic coordinates
- `vars::Matrix{<:AbstractFloat}`: an array of environmental coordinates
"""
struct MikrubiField{T, U <: Real, V <: AbstractFloat}
	ctids::Vector{T}
	locs::VecOrMat{U}
	vars::Matrix{V}
	npixel::Int
	mcounty::Int
	ids::Vector{T}
	starts::Dict{T, Int}
	stops::Dict{T, Int}
	dvar::Int
	function MikrubiField(ctids::AbstractVector{T}, locs::AbstractVecOrMat{U}, 
			vars::AbstractMatrix{V}) where {T, U <: Real, V <: AbstractFloat}
		ctids, locs, vars = copy(ctids), copy(locs), copy(vars)
		npixel = length(ctids)
		npixel == size(locs, 1) == size(vars, 1) ||
			throw(DimensionMismatch(tw"Arguments `ctids`, `locs`, and `vars` 
				should always have the same number of rows!"))
		if ! issorted(ctids)
			perm = sortperm(ctids)
			ctids, locs, vars = ctids[perm], locs[perm, :], vars[perm, :]
		end
		ids = unique(ctids)
		mcounty = length(ids)
		size(vars, 2) >= 5 && 
			@warn tw"It is strongly recommended that no more than four 
				principal components are used for Mikrubi, or parameter space 
				would be highly ill-conditioned!"
		starts = Dict(reverse(ctids) .=> npixel:-1:1)
		stops  = Dict(        ctids  .=> 1:+1:npixel)
		dvar = size(vars, 2)
		new{T, U, V}(ctids, locs, vars, 
			npixel, mcounty, ids, starts, stops, dvar)
	end
end
MikrubiField(ctids, locs, vars) =
	MikrubiField(ctids, colmatrix(locs), colmatrix(float.(vars)))

function Base.show(io::IO, field::MikrubiField)
	print(io, textwrap("Mikrubi Field: 
		geo_dim = $(size(field.locs, 2)),
		env_dim = $(field.dvar),
		$(field.npixel) pixels,
		and $(field.mcounty) counties"))
end

"""
	writefield(path::AbstractString, field::MikrubiField) :: Nothing

Write a Mikrubi field to file at `path`.
"""
function writefield(path::AbstractString, field::MikrubiField)
	headerstring = "I" * "L" ^ size(field.locs, 2) * "V" ^ field.dvar
	header = string.(hcat(headerstring...))
	body = hcat(Array{Any}(field.ctids), field.locs, field.vars)
	writedlm(path, vcat(header, body))
end

"""
	readfield(path::AbstractString) :: MikrubiField

Read a Mikrubi field from file at `path`.
"""
function readfield(path::AbstractString)
	body, header = readdlm(path, Any, header=true)
	heads = header[:]
	Set(heads) == Set(["I", "L", "V"]) && findlast(heads .== "I") == 1 &&
		findfirst(heads .== "V") - findlast(heads .== "L") == 1 ||
			error("The file at `path` is not a well-formatted file!")
	ctids = [body[:, heads .== "I"]...]
	locs = Real.(body[:, heads .== "L"])
	vars = AbstractFloat.(body[:, heads .== "V"])
	MikrubiField(ctids, locs, vars)
end

"""
	MikrubiModel{V <: AbstractFloat}

	MikrubiModel(dvar::Int, params::Vector{<:AbstractFloat})

Construct a Mikrubi Model from a dimensionality `dvar` and a parameter vector
`params`. The equation must hold for `dvar2dparam(dvar) == length(params)`.

Mikrubi Models can be obtained from the function `fit`, and can be used in the
function `predict`.
"""
struct MikrubiModel{V <: AbstractFloat}
	dvar::Int
	params::Vector{V}
	function MikrubiModel(dvar::Int, 
			params::Vector{V}) where {V <: AbstractFloat}
		dvar2dparam(dvar) != length(params) &&
			throw(DimensionMismatch(tw"Length of `params` and value of `dvar` 
				are incompatible! The following equation must hold:
				`(dvar+1) * (dvar+2) / 2 == length(params)`."))
		new{V}(dvar, params)
	end
end

"""
	writemodel(path::AbstractString, model::MikrubiModel) :: Nothing

Write a Mikrubi model to file at `path`.
"""
writemodel(path::AbstractString, model::MikrubiModel) = 
	writedlm(path, vcat(Any[model.dvar], model.params))

"""
	readmodel(path::AbstractString) :: MikrubiModel

Read a Mikrubi model from file at `path`.
"""
function readmodel(path::AbstractString)
	vector = readdlm(path, header=false)[:]
	isinteger(vector[1]) &&
		length(vector)-1 == dvar2dparam(Int(vector[1])) ||
			error("The file at `path` is not a well-formatted file!")
	MikrubiModel(Int(vector[1]), vector[begin+1:end])
end

"""
	writelist(path::AbstractString, list::AbstractVector) :: Nothing

Write any list or vector to file at `path`.
"""
writelist(path::AbstractString, list::AbstractVector) = writedlm(path, list)

"""
	readlist(path::AbstractString) :: Vector

Read any list of vector from file at `path`.
"""
readlist(path::AbstractString) = [readdlm(path, Any, header=false)...]

"""
	decomparams(p::AbstractVector, d::Int) :: Tuple{Matrix, Vector, Any}
	decomparams(model::MikrubiModel) :: Tuple{Matrix, Vector, Any}

Return parameter decomposition `At`, `b`, `c`, where 
- `At` is a lower triangular matrix of size `(d, d)`,
- `b` is a column vector of size `d`, and
- `c` is a scalar.

WARNING: The vector `p` must have length `dvar2dparam(d)`.

# Example
```julia
julia> decomparams(collect(1:10), 3)
([1 0 0; 2 3 0; 4 5 6], [7, 8, 9], 10)
```
"""
function decomparams(p::AbstractVector, d::Int)
	At = zeros(eltype(p), d, d)
	i = j = k = 1
	while i <= d
		At[i, j] = p[k]
		i == j ? i += j = 1 : j += 1
		k += 1
	end
	At, p[k:k+d-1], p[k+d]
end
decomparams(model::MikrubiModel) = decomparams(model.params, model.dvar)

"""
	loglike(field::MikrubiField, params::AbstractVector) 
		:: Vector{<:AbstractFloat}

Compute the log-likelihood of being absent in each pixel given `field` and
`params`.
"""
function loglike(field::MikrubiField, params::AbstractVector)
	At, b, c = decomparams(params, field.dvar)
	pv = field.vars
	loglogistic.(sum((pv * At) .^ 2, dims=2)[:] .+ pv * b .+ c)
end

"""
	findnearest(loc::AbstractVecOrMat{<:Real}, field::MikrubiField) :: Int

Return the row index in `field.locs` which is the nearest to the given 
coordinates.
"""
function findnearest(loc::AbstractVecOrMat{<:Real}, field::MikrubiField)
	loc = loc[:]'
	length(loc) == size(field.locs, 2) ||
		error(textwrap("Dimensionality of `loc` ($(length(loc))) is 
			incompatible with the geographic dimensionality of `field` 
			($(size(field.locs, 2)))!"))
	eucldist2 = sum((loc .- field.locs) .^ 2, dims=2)[:]
	findmin(eucldist2)[2]
end

"""
	findnearests(loc::AbstractVector{<:AbstractVecOrMat}, field::MikrubiField) 
		:: Vector{Int}
	findnearests(loc::AbstractMatrix{<:Real}, field::MikrubiField) 
		:: Vector{Int}

Return the row indices in `field.locs` which are the nearest to each of the 
given coordinates. Duplicate results are reduced to one.
"""
findnearests(loc::AbstractVector{<:AbstractVecOrMat}, field::MikrubiField) =
	unique(findnearest.(loc, [field]))
findnearests(loc::AbstractMatrix, field::MikrubiField) =
	unique(findnearest.(eachrow(float.(loc)), [field]))

"""
	energy(field::MikrubiField, counties, params::AbstractVector)
		:: AbstractFloat
	energy(vars::AbstractMatrix, params::AbstractVector) :: AbstractFloat

Compute the opposite log-likelihood that the occupied counties or occupied 
coordinates are sampled. 
The result is taken opposite sign for optimization, and therefore the function 
is called `energy`.
"""
function energy(field::MikrubiField, counties, params::AbstractVector)
	e = loglike(field, params)
	for o = counties
		start = field.starts[o]
		stop = field.stops[o]
		subsum = sum(e[start:stop])
		e[start:stop] .= 0
		e[start] = log(1 - exp(subsum))
	end
	- sum(e)
end
function energy(vars::AbstractMatrix, params::AbstractVector)
	At, b, c = decomparams(params, size(vars, 2))
	prallp = 1 .- logistic.(sum((vars * At) .^ 2, dims=2)[:] .+ vars * b .+ c)
	- sum(log.(prallp))
end

"""
	fit(field::MikrubiField, counties, coords=zeros(0, 0); 
		optresult=[], iterations=3_000_000, kwargs...) :: MikrubiModel

Numerically find the Mikrubi model maximizing the likelihood that the occupied 
counties as well as the occupied coordinates are sampled in the given Mikrubi 
field. The optimization result is stored in the container `optresult` for 
debugging.
"""
function fit(field::MikrubiField, counties, coords=zeros(0, 0); 
		optresult=[], iterations=39000, kwargs...)
	valcounties = intersect(counties, field.ids)
	indcoords = findnearests(coords, field)
	isempty(valcounties) && isempty(indcoords) &&
		error("No meaningful occupied counties or coordinates!")
	cdvars = field.vars[indcoords, :]
	fun(params) = energy(field, valcounties, params) + energy(cdvars, params)
	# zeroes = zeros(eltype(field.vars), dvar2dparam(field.dvar))
	# TODO: fix the bug in Optim.jl involving NaN32
	zeroes = zeros(Float64, dvar2dparam(field.dvar))
	@info "Now minimizing the opposite likelihood function..."
	result = optimize(fun, zeroes, NelderMead(), Options(
		iterations=iterations, show_trace=true, show_every=500; kwargs...))
	result.iteration_converged && 
		@warn textwrap("The optimizing routine has reached the maximum 
			iteration count (`iterations = $iterations`), and thus the 
			maximizer may be unreliable. Please try to enlarge the parameter.")
	push!(optresult, result)
	@info "Maximized log-likeliness: $(-result.minimum)"
	# MikrubiModel(field.dvar, result.minimizer)
	MikrubiModel(field.dvar, eltype(field.vars).(result.minimizer))
end

"""
	predict(matrix::AbstractMatrix, model::MikrubiModel) :: Vector
	predict(layers::RasterStack, model::MikrubiModel) :: Raster
	predict(field::MikrubiField, model::MikrubiModel) 
		:: Dict{<:Any, <:Vector{<:Tuple{Vector{<:Real}, AbstractFloat}}}

Predict the probability of presence according to processed climatic factors 
(`matrix` / `layers`) or on the Mikrubi `field`.
"""
function predict(matrix::AbstractMatrix, model::MikrubiModel)
	size(matrix, 2) == model.dvar ||
		throw(DimensionMismatch("The number of columns of the matrix 
			($(size(matrix, 2))) is different from the dimensionality of the 
			model ($(model.dvar))!"))
	At, b, c = decomparams(model)
	1 .- logistic.(sum((matrix * At) .^ 2, dims=2)[:] .+ matrix * b .+ c)
end
function predict(layers::RasterStack, model::MikrubiModel)
	matrix, idx = extractlayers(layers)
	layer = makelayer(predict(matrix, model), idx, first(layers))
	return rebuild(layer; name="prob")
end
predict(field::MikrubiField, model::MikrubiModel) =
	Dict(id => predictcounty(field, model, id) for id = field.ids)

"""
	predictcounty(field::MikrubiField, model::MikrubiModel, county) 
		:: Vector{<:Tuple{Vector{<:Real}, AbstractFloat}}

Return the geographic coordinates of pixels in the county sorted by the
likeliness of being occupied.
"""
function predictcounty(field::MikrubiField, model::MikrubiModel, county)
	county in field.ids ||
		error("The Mikrubi field has no such county!")
	idx = field.starts[county] : field.stops[county]
	probs = predict(field.vars[idx, :], model)
	locs = collect.(eachrow(field.locs[idx, :]))
	sort!(tuple.(locs, probs), rev=true, by=last)
end

"""
	probpixels(field::MikrubiField, model::MikrubiModel) 
		:: Vector{<:AbstractFloat}

Compute the probability for every pixel to be occupied in the `field`. 
"""
probpixels(field::MikrubiField, model::MikrubiModel) = 
	predict(field.vars, model)

"""
	probcounties(field::MikrubiField{T, U, V}, model::MikrubiModel{V}) 
		:: Dict{<:Any, <:AbstractFloat}

Compute the probability for every county to be occupied in the `field`. 
"""
function probcounties(field::MikrubiField{T, U, V}, 
		model::MikrubiModel{V}) where {T, U <: Real, V <: AbstractFloat}
	ll = loglike(field, model.params)
	logmprcounties = Dict{T, V}()
	for i = 1:field.npixel
		id = field.ctids[i]
		logmprcounties[id] = ll[i] + get(logmprcounties, id, 0)
	end
	prcounties = Dict{T, V}()
	for id = field.ids
		prcounties[id] = 1 - exp(logmprcounties[id])
	end
	prcounties
end

"""
	samplecounties(field::MikrubiField, model::MikrubiModel) :: Vector{Any}

Sample counties according to their probability of being occupied.
"""
function samplecounties(field::MikrubiField, model::MikrubiModel)
	prcounties = probcounties(field, model)
	@inline bernoulli(p) = rand() <= p
	sort!([k for k = keys(prcounties) if bernoulli(prcounties[k])])
end

"""
	loglipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false) 
		:: AbstractFloat

Calculate the (logarithmic) maximum gradient (in norm) of the probability of 
presence over the `field`. When `wholespace=false` (default), the maximum is 
taken among the points contained in `field`; otherwise it is taken around the 
whole space.
"""
function loglipschitz(model::MikrubiModel, 
		field::MikrubiField; wholespace=false)
	At, b, c = decomparams(model)
	biM = 2 * At * At'
	q(z) = sum((z*At) .^ 2) + sum(z * b) + c
	llpll(qz) = loglogistic(qz) + loglogistic(-qz)
	loglip(z) = llpll(q(z')) + log(sum((biM * z + b).^2))/2
	vars = field.vars
	m, id = findmax(loglip.(eachrow(vars)))
	if wholespace == false
		return m
	else
		return maximum(maximize(loglip, vars[id, :], 
			model.dvar == 1 ? Newton() : NelderMead()))
	end
end

"""
	lipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false) 
		:: AbstractFloat

Calculate the maximum gradient (in norm) of the probability of presence over 
the `field`. When `wholespace=false` (default), the maximum is taken among the 
points contained in `field`; otherwise it is taken around the whole space.
"""
lipschitz(model::MikrubiModel, field::MikrubiField; wholespace=false) = 
	exp(loglipschitz(model, field; wholespace=wholespace))
