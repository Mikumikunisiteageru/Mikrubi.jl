# docs/make.jl

using Mikrubi

using Documenter

makedocs(
	sitename = "A Manual for Mikrubi",
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
