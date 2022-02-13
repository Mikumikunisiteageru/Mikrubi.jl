# example_linear.jl

using Mikrubi

function county_list(R, hw=240)
	sl = 1
	sr = 1 + R
	pixels = -hw+0.5 : hw-0.5
	junctions = vcat(-hw, reverse(-sl:-sl:1-hw), 0:sr:hw-1+sr) .+ 0.5
	cumsum(in.(pixels, (junctions,))), collect(pixels)
end

R = 30
ids, vars = county_list(R)
params = [0.02, 0, 1]
asym = MikrubiField(ids, vars, vars)
model = MikrubiModel(asym, params)

for i = 1:10
	println(sample(model))
end
