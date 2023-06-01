# test/runtests.jl

using Mikrubi
using Aqua
using SafeTestsets

# Aqua.test_ambiguities([Mikrubi, Base, Core])
Aqua.test_unbound_args(Mikrubi)
Aqua.test_undefined_exports(Mikrubi)
Aqua.test_piracy(Mikrubi)
Aqua.test_project_extras(Mikrubi)
Aqua.test_stale_deps(Mikrubi)
Aqua.test_deps_compat(Mikrubi)
Aqua.test_project_toml_formatting(Mikrubi)

@safetestset "functions" begin include("functions.jl") end
@safetestset "rasterizing" begin include("rasterizing.jl") end
@safetestset "exdemo" begin include("exdemo.jl") end
@safetestset "exsim" begin include("exsim.jl") end
@safetestset "exalliwalli" begin include("exalliwalli.jl") end
@safetestset "exjui" begin include("exjui.jl") end
@safetestset "recipesbase" begin include("recipesbase.jl") end
@safetestset "pyplot" begin include("pyplot.jl") end
