# examples/prinsepia/jui.jl

# For Mikrubi.jl, three types of input data are generally required: 
	# (1) a map describing the regions (often a shapefile or a GeoJSON), 
	# (2) raster layers of climatic factors, and 
	# (3) a list of regions occupied by the focal species. 
# The district-level map of Nepal is available from GADM, which can be accessed by commands in Julia using package GADM.jl.

import GADM
shppath = GADM.download("NPL") # ISO code of Nepal

# The climatic factors can likewise be downloaded from WorldClim via the packae RasterDataSources.jl.

using RasterDataSources
ENV["RASTERDATASOURCES_PATH"] = "path/for/data/download"
getraster(WorldClim{BioClim}, res="10m") # 10 min resolution
climpath = RasterDataSources.rasterpath(WorldClim{BioClim})

# The two data sets above should have been prapared in paths `shppath` and `climpath` respectively. Now they can be read and handled with Mikrubi.jl.

using Mikrubi
shptable = readshape(shppath, 1) # District-level
layers = readlayers(climpath)
field, ylayers = makefield(layers, shptable)

# The `makefield` function finishes the map rasterization and the dimensionality reduction of environmental space, and then returns `field`, a structure containing the regionalization and environment information of Nepal, as well as `ylayers`, the low-dimensional environmental factors. 

# According to Pendry (2012), the Prinsepia utilis Royle is present in 17 out of 75 districts of Nepal. The occupied districts are transformed into codes using the function `lookup` and afterwards used to build up a model using `fit`. The model is subsequently applied to `ylayers` using `predict`, which yields the predicted distribution of the species.

regions = ["Dadeldhura", "Doti", "Bajhang", "Kalikot", 
	"Mugu", "Jajarkot", "Jumla", "Rolpa", "Dolpa", "Baglung", 
	"Mustang", "Manang", "Gorkha", "Nuwakot", "Rasuwa", 
	"Okhaldhunga", "Solukhumbu"] # 17 occupied districts
regcodes = lookup(shptable, "NAME_3", regions)
	# Transformed into codes by looking up the NAME_3 column
model = fit(field, regcodes)
geodist = predict(ylayers, model)

# The output `geodist` is of raster format. It can be written to disk using `writelayer` for downstream analysis.

writelayer("path/to/output/geodist.tif", geodist)

# And it can also be illustrated via the package PyPlot.jl (note: this requires Python and its package `matplotlib`) and the submodule `Graphics` of Mikrubi.jl.

using .Graphics
using PyPlot
setplot(PyPlot)
showlayer(geodist)
