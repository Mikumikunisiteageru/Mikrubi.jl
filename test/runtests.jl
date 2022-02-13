using Test
using Mikrubi   # importing the package
cty = ["a", "b", "c", "d", "d", "e", "f", "f", "e"]
                # vector containing county identifiers
geo = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                # matrix of geographic coordinates
env = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                # matrix of environmental coordinates
field = MikrubiField(cty, geo, env)
model = fit(field, ["a", "c", "d"])
                # maximized log-likeliness: -1.9107
# model.pr_cell   # 0.6665, 0.6669, 0.6666, 0.9998, 0.4663, ...
# model.pr_county # "a"=>0.6665, "b"=>0.6669, "c"=>0.6666, ...
# sample(model)   # ["a", "b", "c", "d"] in a random trial
@test isapprox(sum(model.pr_cell), 3.4662445f0)
