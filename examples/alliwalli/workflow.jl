# examples/alliwalli/workflow.jl

using Mikrubi

cd(joinpath(pkgdir(Mikrubi), "examples", "alliwalli"))

# I am not authorized to redistribute the shapefile in the package due to copyright issues. So I show here how the original materials are processed and converted, and their derivatives are attached here as an instance. Users may skip this block when re-running the script. Instead, the derivates are read anew from the disk.
if false
	shptable = readshape("path/to/china/counties.shp");
		# Read a shapefile from Gaode containing all 2894 counties of China
	layers = readlayers("path/to/worldclim/layers");
		# Read 19 climatic factor layers from WorldClim
	china, ylayers = makefield(layers, shptable)
	writefield("chinafield.mkuf", china)
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
china = readfield("chinafield.mkuf")
ylayers = readlayers("ylayers")
ctlist = readlist("countylist.txt")

# Train a model for Allium wallichii in China
model = fit(china, ctlist)
# [ Info: Now minimizing the opposite likelihood function...
# Iter     Function value    √(Σ(yᵢ-ȳ)²)/n
# ------   --------------    --------------
     # 0     4.171830e+04     2.450473e+02
 # * time: 0.003999948501586914
   # 500     1.833867e+02     9.325053e-02
 # * time: 2.2709999084472656
  # 1000     1.470935e+02     1.524631e-01
 # * time: 4.2850000858306885
  # 1500     1.388932e+02     4.105145e-02
 # * time: 6.365000009536743
  # 2000     1.273631e+02     2.085092e-02
 # * time: 8.411999940872192
  # 2500     1.266571e+02     9.378015e-05
 # * time: 10.424999952316284
# [ Info: Maximized log-likeliness: -126.65599400745549
# MikrubiModel{Float64}(3, [1.4842288152354197, -1.3603311815698715, -0.38761691866210646, 1.1231074177981228, 1.2090116395112087, -0.1033479618173679, 14.747024521778938, -14.878922083170924, 11.97056752230023, 30.299436373642205])

# Save the model to disk
writemodel("model.mkum", model)

# Predict the geographic distribution of Allium wallichii in China
geodist = predict(ylayers, model)
# 2160×1080×1 Raster{Float64,3} prob with dimensions:
  # X Projected{Float64} LinRange{Float64}(-180.0, 179.833, 2160) ForwardOrdered Regular Intervals crs: WellKnownText,
  # Y Projected{Float64} LinRange{Float64}(89.8333, -90.0, 1080) ReverseOrdered Regular Intervals crs: WellKnownText,
  # Band Categorical{Int64} 1:1 ForwardOrdered
# extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0), Band = (1, 1))
# missingval: Inf
# crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
# values: [:, :, 1]
           # 89.8333  89.6667  89.5  89.3333  89.1667  89.0  …  -89.1667  -89.3333  -89.5  -89.6667  -89.8333  -90.0
 # -180.0    Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -179.833  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -179.667  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -179.5    Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -179.333  Inf      Inf      Inf   Inf      Inf      Inf   …   Inf       Inf       Inf    Inf       Inf       Inf
 # -179.167  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -179.0    Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
 # -178.833  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
    # ⋮                                        ⋮             ⋱    ⋮                                              ⋮
  # 178.5    Inf      Inf      Inf   Inf      Inf      Inf   …   Inf       Inf       Inf    Inf       Inf       Inf
  # 178.667  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 178.833  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 179.0    Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 179.167  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 179.333  Inf      Inf      Inf   Inf      Inf      Inf   …   Inf       Inf       Inf    Inf       Inf       Inf
  # 179.5    Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 179.667  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf
  # 179.833  Inf      Inf      Inf   Inf      Inf      Inf       Inf       Inf       Inf    Inf       Inf       Inf

# Write the distribution layer to disk
writelayer("geodist.tif", geodist)

# Another application of the function `predict` is to return probabilities organized according to counties
predict(china, model)
# Dict{Int64, Vector{Tuple{Vector{Float64}, Float64}}} with 2893 entries:
  # 1144 => [([113.583, 36.0833], 1.66324e-9), ([113.583, 35.9167], 5.05715e-10), ([113.583, 36.25], 1.03843e-10), ([113. …
  # 2108 => [([107.583, 35.9167], 1.0858e-13), ([107.583, 35.75], 1.05915e-13), ([107.75, 35.75], 4.37428e-14), ([107.75, …
  # 1175 => [([124.917, 46.5833], 0.0), ([124.75, 46.25], 0.0), ([124.917, 46.25], 0.0), ([124.75, 46.0833], 0.0), ([124. …
  # ......

# Calculate the probability of all counties being occupied by the Allium species according to the model
probcounties(china, model)
# Dict{Int64, Float64} with 2893 entries:
  # 2288 => 0.0
  # 1703 => 0.0
  # 1956 => 2.01172e-12
  # 2350 => 1.88738e-15
  # 2841 => 9.8146e-6
  # 2876 => 0.00403172
  # ......

# For the 39-th county (Xiangcheng County, Ganzi Tibetan Autonomous Prefecture, Sichuan Province, China; type locality of Allium xiangchengense), sort its pixels according to probability of being occupied by Allium wallichii
predictcounty(china, model, 39)
# 30-element Vector{Tuple{Vector{Float64}, Float64}}:
 # ([100.08333333333331, 28.583333333333336], 0.014382846487993373)
 # ([99.91666666666663, 28.583333333333336], 0.01280767018326856)
 # ([99.58333333333331, 29.41666666666667], 0.010613283680767083)
 # ([99.41666666666663, 29.41666666666667], 0.010569729644928638)
 # ([99.75, 29.25], 0.01015553778372369)
 # ([99.75, 29.41666666666667], 0.009674145903628695)
 # ......
