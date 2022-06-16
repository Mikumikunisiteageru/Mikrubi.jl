# examples/alliwalli/workflow.jl

cd("examples/alliwalli")

using Mikrubi

# I am not authorized to redistribute the shapefile in the package due to copyright issues. So I show here how the original materials are processed and converted, and their derivatives are attached here as an instance. Users may skip this block when re-running the script. Instead, the derivates are read anew from the disk.
if false
	shptable = readshape("path/to/china/counties.shp");
		# Read a shapefile from Gaode containing all 2894 counties of China
	layers = readlayers("path/to/worldclim/layers");
		# Read 19 climatic factor layers from WorldClim
	china, ylayers = makefield(layers, shptable)
	writefield("chinafield.mikf", china)
		# Save the Mikrubi field to disk
	writelayers("ylayers/ybio*.tif", ylayers)
		# Save the extracted layers to disk
end

# The file `database.tsv` is a result of my georeferencing work for specimens from Chinese Virtual Herbarium (https://www.cvh.ac.cn/) of Allium wallichii. Its first column is taken out and constitutes the list of occupied counties. The list is saved as `countylist.tsv`. When necessary, users may re-run this block. 
if false
	using DelimitedFiles
	db, _ = readdlm("database.tsv", '\t', header=true)
	ctlist = unique(db[:, 1])
	writelist("countylist.txt", ctlist)
end

# Read from disk the Mikrubi field, the extracted layers, and the county list
china = readfield("chinafield.mikf")
ylayers = readlayers("ylayers")
ctlist = readlist("countylist.txt")

# Train a model for Allium wallichii in China
model = fit(china, ctlist)
# [ Info: Maximized log-likeliness: -126.65599400745546
# MikrubiModel{Float64}(3, [1.4842288152354197, -1.3603311815698715, -0.38761691866210646, 1.1231074177981228, 1.2090116395112087, -0.1033479618173679, 14.747024521778938, -14.878922083170924, 11.97056752230023, 30.299436373642205])

# Save the model to disk
writemodel("model.mikm", model)

# Predict the geographic distribution of Allium wallichii in China
geodist = predict(ylayers, model)
# 2160x1080x1 Array{Union{Missing, Float64},3} with AffineMap([0.16666666666666666 0.0; 0.0 -0.16666666666666666], [-180.0, 90.0]) and CRS GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]

# Write the distribution layer to disk
writelayer("geodist.tif", geodist)

# Another application of the function `predict` is to return probabilities organized according to counties
predict(china, model)
# Dict{Int64,Array{Tuple{Array{Float64,1},Float64},1}} with 2893 entries:
  # 2843 => [([102.417, 35.5833], 1.53052e-5), ([102.75, 35.5833], 5.9749e-6), ([102.583, 35.5833], 5.74284e-6), ([102.75 …
  # 1316 => [([117.417, 37.5833], 0.0), ([117.583, 37.5833], 0.0), ([117.75, 37.5833], 0.0), ([117.75, 37.4167], 0.0), ([ …
  # 1333 => [([109.083, 40.4167], 0.0), ([109.25, 40.4167], 0.0), ([109.417, 40.4167], 0.0), ([109.583, 40.4167], 0.0), ( …
  # ......

# Calculate the probability of all counties being occupied by the Allium species according to the model
probcounties(china, model)
# Dict{Int64,Float64} with 2893 entries:
  # 2843 => 3.18779e-5
  # 1316 => 0.0
  # 1333 => 0.0
  # 1671 => 3.48881e-8
  # 1131 => 1.66089e-13
  # 74   => 0.577014
  # ......

# For the 39-th county (Xiangcheng County, Ganzi Tibetan Autonomous Prefecture, Sichuan Province, China; type locality of Allium xiangchengense), sort its pixels according to probability of being occupied by Allium wallichii
predictcounty(china, model, 39)
# 30-element Array{Tuple{Array{Float64,1},Float64},1}:
 # ([100.08333333333331, 28.583333333333336], 0.014382846487993373)
 # ([99.91666666666663, 28.583333333333336], 0.01280767018326856)
 # ([99.58333333333331, 29.41666666666667], 0.010613283680767083)
 # ([99.41666666666663, 29.41666666666667], 0.010569729644928638)
 # ([99.75, 29.25], 0.01015553778372369)
 # ([99.75, 29.41666666666667], 0.009674145903628695)
 # ([99.58333333333331, 28.75], 0.009667323442360876)
 # ......
