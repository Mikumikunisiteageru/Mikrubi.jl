# src/shape.jl

"""
	filterext(dir::AbstractString, extset=nothing) :: Vector{String}

Find all file names in `dir` with extensions in `extset`. When `extset` is set 
to `nothing` (by default), all extensions are acceptable.
"""
function filterext(dir::AbstractString, extset=nothing)
	if isnothing(extset)
		check = isfile
	else
		check = path -> isfile(path) && last(splitext(path)) in extset
	end
	filenames = filter(check, readdir(dir, join=true))
	return filenames
end

"""
	readshape(path::AbstractString, index::Int=-1; 
		extset = [".shp", ".geojson", ".gpkg"]) :: AG.IFeatureLayer

Read a shape file located at `path`. If `path` refers to a file, the file is 
directly read; otherwise, if `path` refers to a directory, a random shape file 
inside is read. 

`extset` describes possible extensions of shape files (see also 
[`readlayers`](@ref)). By setting `extset` to `nothing`, the extension 
filtering is not processed, i.e., all files are regarded as shape files. 
`extset` is indifferent when `path` refers to a file. 

The shape file should contain a dataset. When the dataset consists of multiple 
layers, `index` indicates which data layer should be returned. 
"""
function readshape(path::AbstractString, index::Int=-1; 
		extset = [".shp", ".geojson", ".gpkg"])
	ispath(path) || error("No such file or directory!")
	if isdir(path)
		filenames = filterext(path, extset)
		if length(filenames) == 0
			error("No files with extensions in `extset` in the directory!")
		elseif length(filenames) == 1
			path = first(filenames)
		else
			@warn tw"Multiple files with extensions in `extset` exist in the 
				directory! Now choose a random one to read."
			path = first(filenames)
		end
	end
	dataset = AG.read(path)
	if AG.nlayer(dataset) == 1
		index = 0
	elseif AG.nlayer(dataset) > 1 && index == -1
		display(dataset)
		error(tw"Multiple data layers exist in the dataset! 
			Please designate `index` for one.")
	end
	shptable = AG.getlayer(dataset, index)
	isempty(shptable) && error(tw"No shapes in the dataset!")
	return shptable
end

"""
	goodcolumns(shptable::AG.IFeatureLayer) :: Dict{String, Vector}

Find all properties of features in `shptable` where entries are all unique and 
either integers or strings (types whose `isequal` is well-defined).
"""
function goodcolumns(shptable::AG.IFeatureLayer)
	n = AG.nfield(shptable)
	fields = Dict{String, Vector}()
	feature = first(shptable)
	for i = 0 : AG.nfield(shptable)-1
		list = AG.getfield.(shptable, i)
		if eltype(list) <: Union{Integer, AbstractString} && allunique(list)
			name = AG.getname(AG.getfielddefn(feature, i))
			fields[name] = list
		end
	end
	fields
end

@deprecate goodproperties(shptable) goodcolumns(shptable)

"""
	lookup(shptable::AG.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entry)
	lookup(shptable::AG.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entries::AbstractArray)
	lookup(shptable::AG.IFeatureLayer)

Find row(s) in the shape table whose `column` record(s) equal(s) to `entry` 
or elements of `entries`. When the third argument is an array, results are 
output as an array of the same shape by broadcasting. 
"""
function lookup(shptable::AG.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entry)
	feature = first(shptable)
	i = AG.findfieldindex(feature, column)
	if i == -1 || isnothing(i)
		columns = keys(goodcolumns(shptable))
		error(textwrap("No column in the shapefile named `$column`! 
			Recommended alternatives are $(join(repr.(columns), ", ")). 
			Please select one from them as `column`."))
	end
	shpcol = AG.getfield.(shptable, i)
	index = findall(shpcol .== entry)
	if length(index) == 1
		return index[1]
	elseif length(index) > 1
		@warn textwrap("Multiple entries in the column $column equal to 
			$(repr(entry))! All of these are returned, but please 
			select exactly one from them for succeeding processing.")
		return index
	end # length(index) == 0 afterwards
	if isa(entry, AbstractString) && eltype(shpcol) <: Integer
		@warn tw"Types mismatched: `entry` is a string while the column 
			designated contains integers. Please try again after using `parse` 
			or `tryparse` to convert `entry` to an integer. Nothing is 
			returned here."
		return nothing
	elseif isa(entry, Integer) && eltype(shpcol) <: AbstractString
		@warn tw"Types mismatched: `entry` is an integer while the 
			column designated contains strings. Please try again after using 
			`string` to convert `entry`. Nothing is returned here."
		return nothing
	end
	@warn textwrap("No matched record. Please check the input arguments 
		again. Nothing is returned here.")
	return nothing
end
lookup(shptable::AG.IFeatureLayer, 
		column::Union{AbstractString, Symbol}, entries::AbstractArray) =
	lookup.([shptable], [column], entries)
lookup(shptable::AG.IFeatureLayer) = lookup(shptable, "", 0)
