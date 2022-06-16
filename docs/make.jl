push!(LOAD_PATH, "../src/")

using Mikrubi
using MikrubiGraphics

using Documenter

makedocs(
	sitename = "A Manual for Mikrubi",
	pages = [
		"Home" => "index.md",
		"Manual" => "manual.md",
		"Graphics" => "graphics.md",
		],
	# modules = [Mikrubi, MikrubiGraphics],
)
