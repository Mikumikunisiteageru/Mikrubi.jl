# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

```@meta
CurrentModule = Mikrubi
```

## Reading and writing

### Reading and searching shapefile

Since Mikrubi always reads instead of writes shapefile, only the reading function is offered.

```@docs
readshape
```

The function [`lookup`](@ref) is useful when some attribute (e.g. name or code) of a county is known and the row number of the county in a shapefile is wanted (row number acts as identifiers in the list of occupied counties, see the syntax of [`fit`](@ref)). 

```@docs
lookup
``` 

### Reading and writing list file

List of occupied counties can be prepared explicitly in REPL as a vector or a set. Meanwhile, it is also possible to read from or write to disk such a list, especially when the list is generated outside Julia.

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

It worths mention that when reading layers from a directory, files are sorted according to their names in a manner similar to the sorting order in Windows OS. Please pay extra attention when two parallel series of raster layers are involved for [`makefield`](@ref). 

```@docs
Mikrubi.sortfilenames!
```

All layers read are copied (so that they are no longer read-only and linking to the files on disk) and replaced extreme values by `missing`s. To change this behavior, please refer to the inner function [`Mikrubi.copylayer`](@ref).

```@docs
Mikrubi.copylayer
```

### Reading and writing Mikrubi field

[`MikrubiField`](@ref) is a specially designed type where the environmental information of pixels and their county identifiers are nested. It may be necessary to save (by [`writefield`](@ref)) and load (by [`readfield`](@ref)) a Mikrubi field especially when it is used on multiple species.

```@docs
readfield
writefield
```

### Reading and writing Mikrubi model

[`MikrubiModel`](@ref) is a struct containing transformation parameters. It can be read from and written to disk using respectively [`readmodel`](@ref) and [`writemodel`](@ref). 

```@docs
readmodel
writemodel
```

## Rasterizing a shapefile

Rasterization of one or multiple polygon(s) in a shapefile is implemented in the function [`rasterize`](@ref). The implementation is based on such a fact, that a pixel is touched by a polygon (simple or self-intersecting) when and only when either some edge of the polygon runs through the pixel or the center of the pixel lies inside the polygon. Pixels satisfying the two conditions are collected respectively by internal functions [`Mikrubi.scanline`](@ref) and [`Mikrubi.strokepath`](@ref).

```@docs
rasterize
```

When `rasterize` is applied to an entire shape table (containing multiple counties), its returned value is of an internal type [`Mikrubi.CtPixels`](@ref).

```@docs
Mikrubi.CtPixels
```

#### Internal functions

These internal functions are called directly or indirectly inside [`rasterize`](@ref). 

```@docs
Mikrubi.interpolate
Mikrubi.scanline
Mikrubi.trimedge
Mikrubi.strokeedge
Mikrubi.strokepath
Mikrubi.geom2mat
```

## [Processing the raster layers](@id makefield)

In Mikrubi, climatic factors after being read in typically undergo some processing steps together with shapefile inside the function [`makefield`](@ref), which returns a Mikrubi field and a stack of extracted components in raster layers. The two outputs can be used for training and prediction.

Sometimes it is also required to apply a model to another circumstance (different time or different space), in which case another series of parallel climatic factor layers need to be processed in exactly the same way as those used to generate the Mikrubi field (so that their climatic meanings are the same). Such layers need to be put in the third place in the input argument list for [`makefield`](@ref), and those parallelly extracted components are returned in the third place in output as well.

```@docs
makefield
```

#### Internal functions

These internal functions are called directly or indirectly inside [`makefield`](@ref). Customization is possible when necessary, for example, when the affine transformation during principal component analysis is explicitly wanted.

```@docs
Mikrubi.masklayers!
Mikrubi.extractlayers
Mikrubi.emptylayer
Mikrubi.makelayer
Mikrubi.makelayers
Mikrubi.dftraverse!
Mikrubi.selectvars
Mikrubi.princompvars
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
```
