module Mikrubi

export MikrubiField, logistic, loglogistic, fit, absence, presence

import Optim

Array{T, 2}(a::Array{T, 1}) where {T <: Real} = repeat(a, 1, 1)

struct MikrubiField{T, U <: Real, V <: Real}
	pixel_ids::Array{T, 1}
	pixel_locs::Array{U}
	pixel_vars::Array{V, 2}
	n::Int
	starts::Dict{T, Int}
	stops::Dict{T, Int}
	d::Int
	function MikrubiField(pixel_ids::Array{T}, pixel_locs::Array{U}, 
					  pixel_vars::Array{V}) where {T, U <: Real, V <: Real}
		n = length(pixel_ids)
		n == size(pixel_locs, 1) == size(pixel_vars, 1) ||
			error("Numbers of rows are inconsistent!")
		if ! issorted(pixel_ids)
			perm = sortperm(pixel_ids)
			pixel_ids = pixel_ids[perm]
			pixel_locs = pixel_locs[perm, :]
			pixel_vars = pixel_vars[perm, :]
		end
		starts = Dict(reverse(pixel_ids) .=> n:-1:1)
		stops = Dict(pixel_ids .=> 1:n)
		d = size(pixel_vars, 2)
		new{T, U, V}(pixel_ids, pixel_locs, pixel_vars, n, starts, stops, d)
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
	params = zeros(((field.d+1) * (field.d+2)) >> 1)
	occupieds_ = intersect(occupieds, field.pixel_ids)
	isempty(occupieds_) && error("No meaningful occupied units!")
	fun(params) = energy(field, occupieds_, params)
	result = Optim.optimize(fun, params; kwargs...)
	result.ls_success || println("Warning: Not converged yet!")
	result.minimizer
end

absence(field::MikrubiField, params::Array) = exp.(loglike(field, params))
presence(field::MikrubiField, params::Array) = 1. .- absence(field, params)

end # module
