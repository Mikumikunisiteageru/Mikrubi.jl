# src/Mikrubi.jl

module Mikrubi
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Mikrubi

using Rasters
using RecipesBase

import ArchGDAL; const AG = ArchGDAL
import GeoInterface; const GI = GeoInterface
import DimensionalData; const DD = DimensionalData

import Logistics: Logistic, logistic, loglogistic
import Printf: @sprintf
import DelimitedFiles: readdlm, writedlm
import Statistics: mean, std, cor
import MultivariateStats # fit, projection, PCA
import Optim: optimize, maximize, maximum, NelderMead, Newton, Options
import StatsBase: tiedrank
import ColorTypes: RGB

export logistic, loglogistic
export readshape, goodcolumns, lookup, rasterize
export readlayers, writelayer, writelayers, makefield
export MikrubiField, readfield, writefield
export MikrubiModel, readmodel, writemodel
export readlist, writelist, fit
export predict, predictcounty, probcounties, samplecounties
export lipschitz
export setplot, showlayer, showfield, showctpixels, showshptable

export Graphics # deprecated

const MAXPCADIM = 4

"""
	textwrap(str::AbstractString) :: String

Gobble all linefeeds ("`\\n`") inside `str` and replaces them with spaces
("` `"), so long strings can be wrapped to multiple lines in the codes, like 
the Python package "textwrap". See also [`tw`](@ref).
"""
textwrap(str::AbstractString) = replace(str, r" *\n[ \t]*" => " ")

"""
	@tw_str :: String

Macro version of [`textwrap`](@ref), without interpolation and unescaping.
"""
macro tw_str(str)
	textwrap(str)
end

include("rasterize.jl")
include("shape.jl")
include("layer.jl")
include("core.jl")
include("recipesbase.jl")
include("pyplot.jl")

end # module
