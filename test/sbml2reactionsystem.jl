DATA = joinpath(@__DIR__, "data")
SBML_FILE = joinpath(DATA, "case01.xml")
SYSDICT = sbml_to_sysinfo(SBML_FILE)

# test build_par_map
trueparmap = [Num(Variable{Float64}(:k1)) => 0.8]
parmap = build_par_map(SYSDICT["listOfParameters"])
#=println(typeof(trueparmap[1].first))
println(typeof(parmap[1].first))
println(isequal(trueparmap[1].first, parmap[1].first))=#
@test isequal(parmap, trueparmap)

# test build_comp_map
truecompmap = [Num(Variable{Float64}(:comp1)) => 1.0]
compmap = build_comp_map(SYSDICT["listOfCompartments"])
@test isequal(compmap, truecompmap)

# test build_spec_map
truespecmap = [Num(Variable{Float64}(:A)) => (1.0, Num(Variable{Float64}(:comp1))),
               Num(Variable{Float64}(:B)) => (0.0, Num(Variable{Float64}(:comp1)))]
specmap = build_spec_map(SYSDICT["listOfSpecies"])
@test isequal(specmap, truespecmap)

# test promote_local_parameters
doc = promotelocalparameters(readxml(SBML_FILE))
localpar = lastelement(lastelement(lastelement(lastelement(lastelement(firstelement(doc.root))))))
println(nodename(localpar))
println(localpar["id"])
@test isequal(localpar["id"], "ConvA_k2")
