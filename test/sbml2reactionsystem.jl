DATA = joinpath(@__DIR__, "data")
SBML_FILE = joinpath(DATA, "case01.xml")
SYSDICT = sbml_to_sysinfo(SBML_FILE)

A = Num(Variable{Float64}(:A))
B = Num(Variable{Float64}(:B))
k1 = Num(Variable{Float64}(:k1))
k2 = Num(Variable{Float64}(:k2))
comp1 = Num(Variable{Float64}(:comp1))

# test build_par_map
trueparmap = [k1 => 0.8]
parmap = build_par_map(SYSDICT["listOfParameters"])
@test isequal(parmap, trueparmap)

# test build_comp_map
truecompmap = [comp1 => 1.0]
compmap = build_comp_map(SYSDICT["listOfCompartments"])
@test isequal(compmap, truecompmap)

# test build_spec_map
truespecmap = [A => (1.0, comp1),
               B => (0.0, comp1)]
specmap = build_spec_map(SYSDICT["listOfSpecies"])
@test isequal(specmap, truespecmap)

# test promote_local_parameters
doc = promotelocalparameters(readxml(SBML_FILE))
localpar = lastelement(lastelement(lastelement(lastelement(lastelement(firstelement(doc.root))))))
println(nodename(localpar))
println(localpar["id"])
@test isequal(localpar["id"], "ConvA_k2")

# test getkineticlaw
reaction = firstelement(lastelement(firstelement(readxml(SBML_FILE).root)))
kineticlaw = getkineticlaw(reaction)
@test isequal(nodename(kineticlaw),"kineticLaw")

# test getmath
kineticlaw = lastelement(reaction)
math = getmath(kineticlaw)
@test isequal(nodename(math),"math")

# test _getlistofspeciesreference
reactants, r_stoich = SBML._getlistofspeciesreference(reaction,"Reactant")
products, p_stoich = SBML._getlistofspeciesreference(reaction,"Product")
@test isequal(reactants, [A])
@test isequal(r_stoich, [1])
@test isequal(products, [B])
@test isequal(p_stoich, [1])

# test build_reactions()
reactions = build_reactions([reaction])
kineticlaw = comp1*(A*k1 - (B*k2))
truereactions = [(kineticlaw,[A],[B],[1],[1])]
# @test isequal(reactions[1][1],truereactions[1][1])

# test _process_doc
model = SBML._process_doc(readxml(SBML_FILE))
@test isequal(model.parameters, trueparmap)
@test isequal(model.compartments, truecompmap)
@test isequal(model.species, truespecmap)
@test isequal(model.reactions, truereactions)




#=DATA = "./test/data"
SBML_FILE = joinpath(DATA, "case01.xml")
DOC = readxml(SBML_FILE)
node = DOC.root
ns = namespace(node)
list_name = "listOfParameters"
list_name = "listOfReactions/x:reaction/x:kineticLaw/x:listOfParameters"
list_name = "listOfReactions/x:reaction/x:kineticLaw/x:listOfLocalParameters"
str = "/x:sbml/x:model/x:$(list_name)"
findall(str, node, ["x" => ns])=#