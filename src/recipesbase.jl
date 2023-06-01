# src/recipesbase.jl

pseudorand(x, y) = x % y / y

function hashcolor(x; salt=20, gmin=0.5)
	r0, g0, b0 = pseudorand.(hash((salt, x)), [39, 139, 239])
	r, g, b = @. 1 - (1 - (r0, g0, b0)) * (1 - gmin)
end

function dyecolors(field::MikrubiField; func=tiedrank)
	u(a) = (a .- minimum(a)) ./ (maximum(a) - minimum(a))
	v(a) = u(func(a))
	r = field.dvar >= 1 ? v(field.vars[:, 1]) : fill(0.5, field.npixel)
	g = field.dvar >= 2 ? v(field.vars[:, 2]) : fill(0.5, field.npixel)
	b = field.dvar >= 3 ? v(field.vars[:, 3]) : fill(0.5, field.npixel)
	colores = RGB.(r, g, b)
end

dye(ctpixels::CtPixels; salt=20) = 
	dye(ctpixels, ctpixels.indices; salt=salt)

function dye(ctpixels::CtPixels, grid::Raster; salt=20)
	n, m, _ = size(grid)
	ci = CartesianIndices((n, m))
	imaget = fill(RGB(1.0, 1.0, 1.0), n, m)
	for (ct, pixel) = ctpixels.list
		c = imaget[pixel]
		rt, gt, bt = hashcolor(ct, salt=salt)
		imaget[pixel] = RGB(c.r * rt, c.g * gt, c.b * bt)
	end
	Raster(imaget, dims=dims(grid)[1:2])
end

function dye(field::MikrubiField, grid::Raster; func=tiedrank)
	colores = dyecolors(field; func=func)
	size(field.locs, 2) >= 2 || error("`field.locs` not enough!")
	n, m, _ = size(grid)
	ci = CartesianIndices((n, m))
	imaget = fill(RGB(1.0, 1.0, 1.0), n, m)
	for k = 1:field.npixel
		i, j = xy2ij(grid, field.locs[k, 1:2]...)
		imaget[i, j] = colores[k]
	end
	Raster(imaget, dims=dims(grid)[1:2])
end

function dimends(dim::DD.Dimension)
	dimct = DD.maybeshiftlocus(DD.Center(), dim)
	dimsh = step(dimct) / 2
	return (first(dimct) - dimsh, last(dimct) + dimsh)
end

function dimseq(dim::DD.Dimension)
	dimct = DD.maybeshiftlocus(DD.Center(), dim)
	dimsh = step(dimct) / 2
	ticks = collect(dimct)
	return vcat(ticks .- dimsh, ticks[end] + dimsh)
end

iswhite(rgb::RGB) = isone(rgb)

function xylim(image::AbstractMatrix{<:RGB}, grid::Raster)
	x, y, _... = DD.dims(grid)
	xg0, xg1 = extrema(dimends(x))
	yg0, yg1 = extrema(dimends(y))
	colored = .!iswhite.(image)
	xcolored = any(colored, dims=1)[:]
	ycolored = any(colored, dims=2)[:]
	xl0, xl1 = minmax(dimseq(x)[findfirst(xcolored)], 
					  dimseq(x)[findlast(xcolored)+1])
	yl0, yl1 = minmax(dimseq(y)[findfirst(ycolored)], 
					  dimseq(y)[findlast(ycolored)+1])
	wider = [1.05 -0.05; -0.05 1.05]
	return [xg0, xg1], [yg0, yg1], wider * [xl0, xl1], wider * [yl0, yl1]
end

function xy2ij(layer, x, y)
	DD.dims2indices.(dims(layer)[1:2], [X(Near(x)), Y(Near(y))])
end

struct MPlot end

@recipe f(shptable::AG.IFeatureLayer) = AG.getgeom.(shptable)

@recipe f(ctpixels::CtPixels; salt=20) = MPlot, dye(ctpixels, salt=salt)

@recipe f(grid::Raster, field::MikrubiField; func=tiedrank) = 
	MPlot, dye(field, grid; func=func)

@recipe function f(::Type{MPlot}, mosaic)
	image = Matrix(mosaic')
	xg, yg, xl, yl = xylim(image, mosaic)
	# :xguide --> "X"
	# :yguide --> "Y"
	:xlims --> xl
	:ylims --> yl
	:yflip --> false
	:aspect_ratio --> :auto
	xg, yg, reverse(image, dims=1)
end
