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
                # maximized log-likeliness: -4.3886
model.pr_cell   # 0.6281, 0.4981, 0.3643, 0.4634, 0.3337, ...
model.pr_county # "a"=>0.6281, "b"=>0.4981, "c"=>0.3643, ...
sample(model)   # ["a", "d", "f"] in a random trial
@test isapprox(sum(model.pr_cell), 3.137251f0)
@test isapprox(sum(model.params), -0.7130561f0)
