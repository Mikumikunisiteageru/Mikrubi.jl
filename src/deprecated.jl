# src/deprecated.jl

export Graphics, setplot

module Graphics end

function setplot(x)
	@info "`setplot(PyPlot)` is no longer required for plotting."
	return
end
