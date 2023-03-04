# Mikrubi.jl

```@docs
Mikrubi
```

## A workflow example

We take *Allium wallichii* in China as an example. For more details, please check `examples/alliwalli/workflow.jl`; for the graphic representation of the variables, please visit [Graphics](@ref)).

```julia
julia> using Mikrubi

julia> shptable = readshape(shppath)
Layer: counties
  Geometry 0 (): [wkbPolygon], POLYGON ((99.994226 ...), ...
     Field 0 (id): [OFTInteger64], 2367, 2368, 2369, 2370, 2371, 2372, 2373, ...
     Field 1 (provinceid): [OFTInteger64], 53, 53, 53, 53, 53, 53, 53, 53, ...
     Field 2 (cityid): [OFTInteger64], 5305, 5305, 5305, 5305, 5305, 5323, ...
     Field 3 (cocode): [OFTString], 530524, 530523, 530502, 530521, 530522, ...
     Field 4 (coshname): [OFTString], 昌宁县, 龙陵县, 隆阳区, 施甸县, 腾冲县, 楚雄市, 大姚县, ...
...
 Number of Fields: 14

julia> layers = readlayers(climpath)
[ Info: 19 files "wc2.0_bio_10m_*.tif" recognized in the directory, where * = 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19.
RasterStack with dimensions:
  X Projected{Float64} LinRange{Float64}(-180.0, 179.833, 2160) ForwardOrdered Regular Intervals crs: WellKnownText,
  Y Projected{Float64} LinRange{Float64}(89.8333, -90.0, 1080) ReverseOrdered Regular Intervals crs: WellKnownText,
  Band Categorical{Int64} 1:1 ForwardOrdered
and 19 layers:
  :wc2.0_bio_10m_01 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_02 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_03 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_04 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_05 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_06 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_07 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_08 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_09 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_10 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_11 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_12 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_13 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_14 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_15 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_16 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_17 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_18 Float64 dims: X, Y, Band (2160×1080×1)
  :wc2.0_bio_10m_19 Float64 dims: X, Y, Band (2160×1080×1)

julia> ctlist = readlist(ctlistpath)
46-element Array{Int64,1}:
  568
  162
  364
  ...
  233
 2768
 2770

julia> field, ylayers = makefield(layers, shptable);

julia> field
Mikrubi Field: geo_dim = 2, env_dim = 3, 62716 pixels, and 2893 counties

julia> ylayers
RasterStack with dimensions:
  X Projected{Float64} LinRange{Float64}(-180.0, 179.833, 2160) ForwardOrdered Regular Intervals crs: WellKnownText,
  Y Projected{Float64} LinRange{Float64}(89.8333, -90.0, 1080) ReverseOrdered Regular Intervals crs: WellKnownText,
  Band Categorical{Int64} 1:1 ForwardOrdered
and 3 layers:
  :pca1 Float64 dims: X, Y, Band (2160×1080×1)
  :pca2 Float64 dims: X, Y, Band (2160×1080×1)
  :pca3 Float64 dims: X, Y, Band (2160×1080×1)

julia> model = fit(field, ctlist)
[ Info: Now minimizing the opposite likelihood function...
Iter     Function value    √(Σ(yᵢ-ȳ)²)/n
------   --------------    --------------
     0     4.171830e+04     2.450473e+02
 * time: 0.01399993896484375
   500     1.833867e+02     9.325053e-02
 * time: 2.998000144958496
  1000     1.470935e+02     1.524631e-01
 * time: 4.9700000286102295
  1500     1.388932e+02     4.105145e-02
 * time: 6.976000070571899
  2000     1.273631e+02     2.085092e-02
 * time: 8.812000036239624
  2500     1.266571e+02     9.378015e-05
 * time: 10.5239999294281
[ Info: Maximized log-likeliness: -126.65599400745549
MikrubiModel{Float64}(3, [1.4842288152354197, -1.3603311815698715, -0.38761691866210646, 1.1231074177981228, 1.2090116395112087, -0.1033479618173679, 14.747024521778938, -14.878922083170924, 11.97056752230023, 30.299436373642205])

julia> geodist = predict(ylayers, model)
2160×1080×1 Raster{Float64,3} prob with dimensions:
  X Projected{Float64} LinRange{Float64}(-180.0, 179.833, 2160) ForwardOrdered Regular Intervals crs: WellKnownText,
  Y Projected{Float64} LinRange{Float64}(89.8333, -90.0, 1080) ReverseOrdered Regular Intervals crs: WellKnownText,
  Band Categorical{Int64} 1:1 ForwardOrdered
extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0), Band = (1, 1))
missingval: -1.7e308
crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
values: [:, :, 1]
           89.8333   89.6667   89.5      89.3333   89.1667   …  -89.3333   -89.5      -89.6667   -89.8333   -90.0
 -180.0    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.833  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.667  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.5    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.333  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308  …   -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.167  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -179.0    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
 -178.833  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
    ⋮                                               ⋮        ⋱                                                ⋮
  178.5    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308  …   -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  178.667  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  178.833  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.0    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.167  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.333  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308  …   -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.5    -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.667  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308
  179.833  -1.7e308  -1.7e308  -1.7e308  -1.7e308  -1.7e308      -1.7e308   -1.7e308   -1.7e308   -1.7e308   -1.7e308

julia> writelayer("path/to/output/geodist.tif", geodist)

```

## Outline of [Manual](@ref)

```@contents
Pages = [
    "manual.md",
]
Depth = 3
```

## Outline of [Graphics](@ref)

```@contents
Pages = [
    "graphics.md",
]
Depth = 4
```

## [Index](@id main-index)

```@index
```
