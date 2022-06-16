# examples/onedimsim/sim.jl

using Mikrubi

# Prepare inputs for skewed field where counties have different sizes
function county_list(R, hw=240)
	sl = 1      # Size of counties on the left size (negative)
	sr = 1 + R  # Size of counties on the right size (positive)
	pixels = -hw+0.5 : hw-0.5
	junctions = vcat(-hw, reverse(-sl:-sl:1-hw), 0:sr:hw-1+sr) .+ 0.5
	cumsum(in.(pixels, (junctions,))), collect(pixels)
end

# Set the parameter R, and generate required arguments
R = 30
ctids, vars = county_list(R)

# Fix the parameters
params = [0.02, 0, 1]

# Construct the field directly
asym = MikrubiField(ctids, vars, vars)

# Construct the model directly
model = MikrubiModel(1, params)

# Sample the counties for ten times
for _ = 1:10
	println(samplecounties(asym, model))
end
