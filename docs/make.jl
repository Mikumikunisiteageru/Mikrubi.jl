# docs/make.jl

using Mikrubi

using Documenter
using PyPlot

makedocs(
	sitename = "Mikrubi.jl",
	pages = [
		"Home" => "index.md",
		"Manual" => "manual.md",
		"Graphics" => "graphics.md",
		],
	modules = [Mikrubi],
)

deploydocs(
    repo = "github.com/Mikumikunisiteageru/Mikrubi.jl.git",
	versions = ["stable" => "v^", "v#.#.#"]
)
