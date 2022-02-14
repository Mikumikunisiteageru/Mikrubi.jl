using Test
using Mikrubi   # importing the package
cty = ["a", "b", "c", "d", "d", "e", "f", "f", "e"]
                # vector containing county identifiers
geo = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                # matrix of geographic coordinates
env = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                # matrix of environmental coordinates
field = MikrubiField(cty, geo, env)
model = fit(field, ["a", "c", "f"])
                # maximized log-likeliness: -4.389
model.pr_pixel  # 0.628, 0.498, 0.364, 0.463, 0.334, 0.222 ...
model.pr_county # "a" => 0.628, "b" => 0.498, "c" => 0.364 ...
sample(model)   # ["a", "d", "f"] in a random trial
@test isapprox(sum(model.pr_pixel), 3.137251f0)
@test isapprox(sum(model.params), -0.7130561f0)
