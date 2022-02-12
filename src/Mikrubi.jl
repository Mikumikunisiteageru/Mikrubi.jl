module Mikrubi

export MikrubiField, logistic, loglogistic, fit, MikrubiModel, sample

import Optim

Matrix{T}(a::Vector{T}) where {T <: Real} = repeat(a, 1, 1)

struct MikrubiField{T, U <: Real, V <: AbstractFloat}
	pixel_ids::Array{T, 1}
	pixel_locs::Array{U}
	pixel_vars::Array{V, 2}
	n::Int
	m::Int
	ids::Vector{T}
	starts::Dict{T, Int}
	stops::Dict{T, Int}
	d::Int
	function MikrubiField(pixel_ids::Array{T}, pixel_locs::Array{U}, 
			pixel_vars::Array{V}) where {T, U <: Real, V <: AbstractFloat}
		n = length(pixel_ids)
		n == size(pixel_locs, 1) == size(pixel_vars, 1) ||
			error("Numbers of rows are inconsistent!")
		if ! issorted(pixel_ids)
			perm = sortperm(pixel_ids)
			pixel_ids = pixel_ids[perm]
			pixel_locs = pixel_locs[perm, :]
			pixel_vars = pixel_vars[perm, :]
		end
		ids = unique(pixel_ids)
		m = length(ids)
		starts = Dict(reverse(pixel_ids) .=> n:-1:1)
		stops = Dict(pixel_ids .=> 1:n)
		d = size(pixel_vars, 2)
		new{T, U, V}(pixel_ids, pixel_locs, pixel_vars, n, m, ids,
			starts, stops, d)
	end
end

logistic(x :: Real) = 1. ./ (1. + exp(-x))
loglogistic(x :: Real) = -log(1. + exp(-x))

function param_decompose(p::Array, d::Int)
	A = zeros(d, d)
	i = j = k = 1
	while i <= d
		A[i, j] = p[k]
		i == j ? i += j = 1 : j += 1
		k += 1
	end
	A, p[k:k+d-1], p[k+d]
end

function loglike(field::MikrubiField, params::Array)
	A, b, c = param_decompose(params, field.d)
	pv = field.pixel_vars
	loglogistic.(sum((pv * A) .^ 2, dims=2)[:] .+ pv * b .+ c)
end

function energy(field::MikrubiField, occupieds::Union{Array, Set}, params::Array)
	e = loglike(field, params)
	for o = occupieds
		start = field.starts[o]
		stop = field.stops[o]
		subsum = sum(e[start:stop])
		e[start:stop] .= 0.
		e[start] = log(1. - exp(subsum))
	end
	- sum(e)
end

function fit(field::MikrubiField, occupieds::Union{Array, Set}; kwargs...)
	occupieds_ = intersect(occupieds, field.ids)
	isempty(occupieds_) && error("No meaningful occupied units!")
	fun(params) = energy(field, occupieds_, params)
	zeroes = zeros(((field.d+1) * (field.d+2)) >> 1)
	result = Optim.optimize(fun, zeroes, iterations=3000000; kwargs...)
	result.iteration_converged && println("Warning: Not converged yet!")
	MikrubiModel(field, result.minimizer)
end

struct MikrubiModel{T, U <: Real, V <: AbstractFloat}
	field::MikrubiField{T, U, V}
	params::Vector{V}
	pr_cell::Vector{V}
	pr_county::Dict{T, V}
	function MikrubiModel(field::MikrubiField{T, U, V}, 
			params::Vector{V}) where {T, U <: Real, V <: AbstractFloat}
		ll = loglike(field, params)
		pr_cell = 1. .- exp.(ll)
		t = Dict{T, V}()
		for i = 1:field.n
			id = field.pixel_ids[i]
			t[id] = ll[i] + get(t, id, 0)
		end
		pr_county = Dict{T, V}()
		for id = field.ids
			pr_county[id] = 1. - exp(t[id])
		end
		new{T, U, V}(field, params, pr_cell, pr_county)
	end
end

function sample(model::MikrubiModel)
	pr_county = model.pr_county
	@inline bernoulli(p) = rand() <= p
	sort!([k for k = keys(pr_county) if bernoulli(pr_county[k])])
end

end # module
