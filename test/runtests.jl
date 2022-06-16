using Test
using Mikrubi                 # importing the package
cty = ["a", "b", "c", "d", "d", "e", "f", "f", "e"]
                               # vector containing county identifiers
geo = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                               # matrix of geographic coordinates
env = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3]
                               # matrix of environmental coordinates
field = MikrubiField(cty, geo, env)
model = fit(field, ["a", "c", "f"])
                               # maximized log-likeliness: -4.389
samplecounties(field, model)   # ["a", "d", "f"] in a random trial
@test isapprox(sum(Mikrubi.probpixels(field, model)), 3.137251f0)
@test isapprox(sum(values(probcounties(field, model))), 2.8943691f0)
