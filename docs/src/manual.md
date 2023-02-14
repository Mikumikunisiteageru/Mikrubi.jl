# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

```@meta
Module = Mikrubi
```

## Reading and writing

### Reading and searching shapefile

Since shape files are always read instead of written in this application, only the reading function [`readshape`](@ref) is provided. 

```@docs
readshape
```

The function [`lookup`](@ref) is useful when some attribute (e.g. name or code) of a county is known and the row number of the county in a shapefile is wanted (row number may act as identifiers in the list of occupied counties, see the syntax of [`fit`](@ref)). 

```@docs
lookup
``` 

#### Internal functions

```@docs
Mikrubi.filterext
Mikrubi.goodcolumns
```

### Reading and writing list file

List of occupied counties can be prepared explicitly in Julia as a vector or a set. Meanwhile, it is also possible to read from or write to disk such a list, especially when the list is generated outside Julia.

```@docs
readlist
writelist
```

### Reading and writing raster layers

Climatic factors are downloaded and stored as raster layers. Mikrubi reads such layers by [`readlayers`](@ref), performs principal component analysis on them and returns the results as layers also. When the output layers need to be kept for future use, they can be written to disk using [`writelayers`](@ref). Moreover, when the predicted distribution of species is organized in raster format, it can be saved likewise using [`writelayer`](@ref).

```@docs
readlayers
writelayer
writelayers
```

#### Internal functions

It is worth mention that when reading layers from a directory, files are sorted according to their names in a manner similar to the sorting order in Windows OS. Please pay extra attention when two parallel raster stacks are fed into [`makefield`](@ref). 

```@docs
Mikrubi.sortfilenames!
Mikrubi.allsame
```

### Reading and writing Mikrubi fields

[`MikrubiField`](@ref) is a specially designed type where the environmental information of pixels and their county identifiers are nested. It may be necessary to save (by [`writefield`](@ref)) and load (by [`readfield`](@ref)) a Mikrubi field especially when it is used on multiple species.

```@docs
readfield
writefield
```

### Reading and writing Mikrubi models

[`MikrubiModel`](@ref) is a struct containing transformation parameters. It can be read from and written to disk using respectively [`readmodel`](@ref) and [`writemodel`](@ref). 

```@docs
readmodel
writemodel
```

## Rasterizing a shapefile

Since v1.3.0, Mikrubi no longer provides its own rasterization routine; the implementation from [Rasters](https://rafaqz.github.io/Rasters.jl/stable/) is applied instead. The function [`rasterize`](@ref) in Mikrubi integrates the rasterization of multiple geometries. The returned value is of an internal type [`Mikrubi.CtPixels`](@ref).

```@docs
rasterize
Mikrubi.CtPixels
```

#### Internal functions

```@docs
Mikrubi.getpixels
Mikrubi.getcounties
Mikrubi.getpixel
Mikrubi.getcounty
Mikrubi.indicate
Mikrubi.register!
Mikrubi.ispoly
```

## [Processing the raster layers](@id makefield)

In Mikrubi, climatic factors after being read in typically undergo some processing steps together with shapefile inside the function [`makefield`](@ref), which returns a Mikrubi field and a stack of extracted components in raster layers. The two outputs can be used for training and prediction.

Sometimes it is also required to apply a model to another circumstance (different time or different space), in which case another series of parallel climatic factor layers need to be processed in exactly the same way as those used to generate the Mikrubi field (so that their climatic meanings are the same). Such layers need to be put in the third place in the input argument list for [`makefield`](@ref), and those parallelly extracted components are returned in the third place in output as well.

```@docs
makefield
```

#### Internal functions

```@docs
Mikrubi.colmatrix
Mikrubi.masklayers!
Mikrubi.extractlayers
Mikrubi.emptylayer!
Mikrubi.emptylayer
Mikrubi.emptylayers
Mikrubi.makelayer
Mikrubi.makelayers
Mikrubi.dftraverse!
Mikrubi.selectvars
Mikrubi.princompvars
Mikrubi.DimLower
Mikrubi.dimpoints
Mikrubi.centercoords
Mikrubi.buildfield
```

## The Mikrubi core

Two specially designed structs are involved in the core of Mikrubi.

### Mikrubi field

[`MikrubiField`](@ref) is a struct containing mainly three types of information of pixels/points, that is, which counties they belong to (`ctids`), their geographic coordinates (`locs`), and their environmental coordinates (`vars`), with also some derived assistant attributes, such as geographic dimensionality (usually `2`) and environmental dimensionality (for example, `3`).

[`MikrubiField`](@ref) can be obtained in three ways: as first output argument of [`makefield`](@ref makefield), read from disk, or constructed directly from the three required attributes (this may be useful for simulation analysis). 

```@docs
MikrubiField
```

### Mikrubi model

[`MikrubiModel`](@ref) contains the environmental dimensionality and the model parameters to define a positive-definite quadratic mapping from environmental space to a real number axis. 

Like [`MikrubiField`](@ref), [`MikrubiField`](@ref) can be obtained in three ways: as output argument of [`fit`](@ref fitfuncm), read from disk, or constructed directly from attributes. An example of obtaining Mikrubi field and Mikrubi model from constructors are available in `examples/onedimsim/sim.jl`.

```@docs
MikrubiModel
```

### [Fitting a Mikrubi model](@id fitfuncm)

When a Mikrubi field as well as occurrence data in county and/or in coordinates are ready, they can be used to train a Mikrubi model by function [`fit`](@ref) (county data at `counties`, required; coordinates at `coords`, optional). Result is output as a [`MikrubiModel`](@ref).

```@docs
fit(field::MikrubiField, counties, coords=zeros(0, 0))
```

### Predicting from a Mikrubi model

A Mikrubi model can be applied by function [`predict`](@ref) to a matrix with its columns corresponding to extracted variables, a stack of extracted layers, or a Mikrubi field. 
- When input argument is a matrix, output argument is a column vector denoting the probability of presence in pixels/points related to rows in the matrix.
- When input argument is a stack of layers, output argument is a single layer denoting the probability of presence.
- When input argument is a Mikrubi field, output argument is a `Dict` which maps every county identifier to probability of presence at pixels inside the county, see also [`predictcounty`](@ref).

```@docs
predict
```

When distribution probability within only one county is concerned, [`predictcounty`](@ref) returns probability of presence at all pixels that constitute the county in descending order. Therefore, the first element represents the most likely occupied pixel of a county.

```@docs
predictcounty
```

It is also possible to obtain the overall probability that every county is occupied by the function [`probcounties`](@ref).

```@docs
probcounties
```

### Sampling counties in a Mikrubi field

For simulation analysis, sometimes it is required to sample a set of counties from a Mikrubi field and a Mikrubi model. [`samplecounties`](@ref) does the trick.

```@docs
samplecounties
```

### Detecting overfitting

Overfitting can be detected with the Lipschitz constant, the (logarithmic) maximum gradient (in norm) of the probability of presence in environmental space.

```@docs
lipschitz
```

### Mathematic functions

Two mathematic functions used in the core are exported for convenience.

```@docs
logistic
loglogistic
```

#### Internal functions

```@docs
Mikrubi.dvar2dparam
Mikrubi.decomparams
Mikrubi.loglike
Mikrubi.energy
Mikrubi.probpixels
Mikrubi.findnearest
Mikrubi.findnearests
Mikrubi.loglipschitz
Mikrubi.textwrap
Mikrubi.@tw_str
```
