# Mikrubi.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://Mikumikunisiteageru.github.io/Mikrubi.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://Mikumikunisiteageru.github.io/Mikrubi.jl/dev)
[![CI](https://github.com/Mikumikunisiteageru/Mikrubi.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Mikumikunisiteageru/Mikrubi.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/Mikumikunisiteageru/Mikrubi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Mikumikunisiteageru/Mikrubi.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Mikrubi Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Mikrubi)](https://pkgs.genieframework.com?packages=Mikrubi)

*Mikrubi: a model for species distributions using region-based records*

Many species occurrence records from specimens and publications are based on regions such as administrative units (thus sometimes called `counties` in the codes). These region-based records are accessible and dependable, and sometimes they are the only available data source; however, few species distribution models accept such data as direct input. In [Yang et al. (2023)](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06283), we present a method named Mikrubi for robust prediction of species distributions from region-based occurrence data. This is the Julia package implementing the algorithms. 

## Installation

Mikrubi.jl currently requires Julia v1.7.0 or higher. This registered package can be installed inside the Julia REPL by typing
```julia
]add Mikrubi
```

## Input data requirements

To estimate the fine-scale distribution of a species using its presence or absence in each region, the package generally requires three types of input data: 

- A map describing the shapes of all regions as polygons. For many countries or regions, such an administrative partition map can be found from [Database of Global Administrative Areas](https://gadm.org/) (accessible via [GADM.jl](https://github.com/JuliaGeo/GADM.jl), see examples/prinsepia/jui.jl). Specifically for China, the correct county-level shapefile is available from [National Platform of Common Geospatial Information Services](https://www.tianditu.gov.cn/) and [Gaode Map Open Platform](https://lbs.amap.com/).

- Raster layers of climatic factors of the same size, shape, and resolution. A commonly used dataset is [WorldClim](https://worldclim.org/data/index.html) (accessible via [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl)). 

- A list of regions occupied by the species. 

## Workflow

A typical workflow of the package resembles the following lines, where `shppath` refers to the path to the map file, `climpath` refers to the directory path to the raster files, and `ctlistpath` refers to the path to the list containing lines of integer identifiers representing the regions.

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

## Citation

An introduction of this package and the model it implements has been published on [Ecography (10.1111/ecog.06283)](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06283).

Thank Michael Krabbe Borregaard @mkborregaard and Rafael Schouten @rafaqz for reviewing the code and the manuscript and giving greatly constructive opinions!

If you apply the package or the model in your research, please cite them via the paper above or as the following after substituting the version:
```
Yang, Y.-C., Zhang, Q. and Chen, Z.-D. 2023. Mikrubi: a model for species distributions using region-based records. – Ecography 2023: e06283 (ver. 1.3.2).
```
The equivalent BibTeX file for citation is available at [CITATION.bib](https://github.com/Mikumikunisiteageru/Mikrubi.jl/blob/master/CITATION.bib). You may also import this file or the following BibTeX code block to your reference management software.
```bibtex
@article{Mikrubi2023,
	author = {Yang, Yu-Chang and Zhang, Qian and Chen, Zhi-Duan},
	title = {Mikrubi: a model for species distributions using region-based records},
	journal = {Ecography},
	year = {2023},
	volume = {2023},
	pages = {e06283},
	doi = {https://doi.org/10.1111/ecog.06283},
	url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.06283},
	eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/ecog.06283},
}
```

#### Patch to the Ecography paper

Due to [an update of GADM.jl on May 16, 2023](https://github.com/JuliaGeo/GADM.jl/commit/f7bebc9c358a9d00540e42e90e47ad7b6ca145bf) which fetches GADM data v4.1 rather than v3.6, the district-level map of Nepal now has index `3` (no longer `1` as in the paper). The corresponding line in the code example from the paper (see also [jui.jl](https://github.com/Mikumikunisiteageru/Mikrubi.jl/blob/master/examples/prinsepia/jui.jl)) should be changed to

```julia
shptable = readshape(shppath, 3) # District-level
```
