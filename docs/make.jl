# docs/make.jl

using Mikrubi
using .Graphics

cd(joinpath(pkgdir(Mikrubi), "docs"))

using Documenter

makedocs(
	sitename = "A Manual for Mikrubi",
	pages = [
		"Home" => "index.md",
		"Manual" => "manual.md",
		"Graphics" => "graphics.md",
		],
	# modules = [Mikrubi, Mikrubi.Graphics],
)
