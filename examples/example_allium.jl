# example_allium.jl

using DelimitedFiles
using Mikrubi

data = readdlm("field_china.tsv")
china = MikrubiField(Int.(data[:, 1]), data[:, 2:3], data[:, 4:6])
allium = [2858, 68, 22, 17, 72, 2859, 16, 2770, 33, 233, 89, 73, 18, 
          2768, 74, 59, 92, 34, 69, 3, 35, 19, 38, 2850, 79, 2765, 2853, 
          2864, 1010, 82, 216, 36, 786, 162, 568, 364, 386, 454, 291, 48, 
          2744, 2740, 2507, 98, 134, 488] # list of occupied counties
model = fit(china, allium) # needs about 15 s, log-likeliness = -126.6560
prop = model.pr_pixel

topo = fill(NaN, 1080, 540)
c2i(x) = floor(Int, x * 6) # coords to index
for i = 1:china.n
	x = c2i(china.pixel_locs[i, 1])
	y = c2i(china.pixel_locs[i, 2])
	topo[x, y] = prop[i]
end
crop = topo[c2i(73+.01):c2i(134-.01), c2i(54-.01):-1:c2i(10+.01)]

using PyPlot
imshow(crop'.^0.4, extent=(73, 134, 10, 54))
gca().set_aspect(sec(pi/6))
set_cmap("CMRmap")
