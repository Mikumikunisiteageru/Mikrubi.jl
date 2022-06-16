# Mikrubi.jl

*A model for species distributions on county-level records.*

Mikrubi is a species distribution model (SDM) based on rough occurrences and/or precise coordinates. 

If you want to model the geographic distribution or the ecological niche of a species on its county-level records (or in any other administrative unit) and maybe also with some coordinates, Mikrubi is be a good option.

## Workflow

Suppose we want to model the fine-scale distribution of a certain species while knowing its occupied counties of a certain country. To start Mikrubi, we need to prepare three requisites:

- A shapefile located at the path `shpfile`, describing shapes of all counties as polygons. For many countries or regions, such an administrative partition shapefile can be found from [Database of Global Administrative Areas](https://gadm.org/). Especially for China, the accepted county-level shapefile is available from [National Platform of Common Geospatial Information Services](https://www.tianditu.gov.cn/) and [Gaode Map Open Platform](https://lbs.amap.com/).

- A directory `layerdir`, containing typically multiple raster layers usually from [WorldClim](https://worldclim.org/data/index.html), of the same size, shape, and resolution.  

- A list of integers explicitly provided or written in file at path `countylist`, each integer of which corresponds to the row number of an occupied county in the `shpfile`. 

Then the workflow can be summarized as the following:

```julia
using Mikrubi
shptable = readshape(shpfile)
layers = readlayers(layerdir)
field, ylayers = makefield(layers, shptable)
ctlist = readlist(countylist)
model = fit(field, ctlist)
geodist = predict(ylayers, model)
writelayer("path/to/output/geodist.tif", geodist)
```

Now we take *Allium wallichii* in China as an example (for the same case, codes with more details are also available in `examples/alliwalli/workflow.jl`; for graphic representation of variables here, please visit [Mikrubi Graphics](@ref)):

```julia
julia> using Mikrubi

julia> shptable = readshape(shpfile)
Shapefile.Table{Union{Missing, Shapefile.Polygon}} with 2894 rows and the following 15 columns:

geometry, id, provinceid, cityid, cocode, coshname, prcode, cishname, cicode, couname, prshname, cifullname, prfullname, countyid, area


julia> layers = readlayers(layerdir)
[ Info: Totally 19 raster files recognized in the directory, in the form of wc2.0_bio_10m_***.tif, where asterisks mean 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19.
19-element Array{GeoArray{Union{Missing, Float64}},1}:
 2160x1080x1 Array{Union{Missing, Float64},3} with AffineMap([0.16666666666666666 0.0; 0.0 -0.16666666666666666], [-180.0, 90.0]) and CRS GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
 ... (18 more omitted)

julia> field, ylayers = makefield(layers, shptable);
[ Info: Rasterizing procedure may be slow. Please wait...
################################################################################
[ Info: Among the 63085 pixels, 369 (0.6%) is/are discarded for lacking values.

julia> field
Mikrubi Field: geo_dim = 2, env_dim = 3, 62716 pixels, and 2893 counties

julia> ylayers
3-element Array{GeoArray{Union{Missing, Float64}},1}:
 2160x1080x1 Array{Union{Missing, Float64},3} with AffineMap([0.16666666666666666 0.0; 0.0 -0.16666666666666666], [-180.0, 90.0]) and CRS GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
 ... (2 more omitted)

julia> ctlist = readlist(countylist)
46-element Array{Int64,1}:
  568
  162
  364
  ...
  233
 2768
 2770

julia> model = fit(field, ctlist)
[ Info: Maximized log-likeliness: -126.65599400745546
MikrubiModel{Float64}(3, [1.4842288152354197, -1.3603311815698715, -0.38761691866210646, 1.1231074177981228, 1.2090116395112087, -0.1033479618173679, 14.747024521778938, -14.878922083170924, 11.97056752230023, 30.299436373642205])

julia> geodist = predict(ylayers, model)
2160x1080x1 Array{Union{Missing, Float64},3} with AffineMap([0.16666666666666666 0.0; 0.0 -0.16666666666666666], [-180.0, 90.0]) and CRS GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]

julia> writelayer("path/to/output/geodist.tif", geodist)
"path/to/output/geodist.tif"

```

## Manual Outline

```@contents
Pages = [
    "manual.md",
]
Depth = 3
```

## Graphics Outline

```@contents
Pages = [
    "graphics.md",
]
Depth = 3
```

## [Index](@id main-index)

```@index
```
