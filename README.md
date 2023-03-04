# Mikrubi.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://Mikumikunisiteageru.github.io/Mikrubi.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Mikumikunisiteageru.github.io/Mikrubi.jl/dev)
[![CI](https://github.com/Mikumikunisiteageru/Mikrubi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/Mikumikunisiteageru/Mikrubi.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/Mikumikunisiteageru/Mikrubi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Mikumikunisiteageru/Mikrubi.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

*Mikrubi: a model for species distributions using region-based records*

Many species occurrence records from specimens and publications are based on regions such as administrative units (thus sometimes called `counties` in the codes). These region-based records are accessible and dependable, and sometimes they are the only available data source; however, few species distribution models accept such data as direct input. In Yang et al. (unpublished), we present a method named Mikrubi for robust prediction of species distributions from region-based occurrence data. This is the Julia package implementing the algorithms. 

## Input data requirements

To estimate the fine-scale distribution of a species using its presence or absence in each region, the package generally requires three types of input data: 

- A map describing the shapes of all regions as polygons. For many countries or regions, such an administrative partition map can be found from [Database of Global Administrative Areas](https://gadm.org/) (accessible via [GADM.jl](https://github.com/JuliaGeo/GADM.jl), see examples/prinsepia/jui.jl). Specifically for China, the correct county-level shapefile is available from [National Platform of Common Geospatial Information Services](https://www.tianditu.gov.cn/) and [Gaode Map Open Platform](https://lbs.amap.com/).

- Raster layers of climatic factors of the same size, shape, and resolution. A commonly used dataset is [WorldClim](https://worldclim.org/data/index.html) (accessible via [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl)). 

- A list of regions occupied by the species. 

## Workflow

A typical workflow of the package resembles the following lines, where `shppath` refers to the path to the map file, `climpath` refers to the directory path to the raster files, and `ctlistpath` refers to the path to the list containing lines of integers representing the regions.

```julia
using Mikrubi
shptable = readshape(shppath)
layers = readlayers(climpath)
ctlist = readlist(ctlistpath)
field, ylayers = makefield(layers, shptable)
model = fit(field, ctlist)
geodist = predict(ylayers, model)
writelayer("path/to/output/geodist.tif", geodist)
```
