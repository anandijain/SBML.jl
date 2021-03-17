using Test, SBML, EzXML, MathML, Symbolics

fn = "data/01193-sbml-l3v1.xml"
d = sbml_to_sysinfo(fn)
@test isa(d, Dict{String,Union{Nothing,Vector{EzXML.Node}}})
@test !all(isnothing, values(d))

fn = "data/BIOMD0000000346_url.xml"
xml = readxml(fn)
d = sbml_to_sysinfo(fn)
td = make_nice(d)
@test isa(td, Dict)

# Namespacing issue: fbc:ChemicalFormula
@test isempty(collect(skipmissing(td["listOfSpecies"][end].second)))

# start of ReactionNetwork bind

# https://github.com/SciML/CellMLToolkit.jl/blob/dev/src/cellml.jl
# * make find_iv work
fn = "data/BIOMD0000000704_url.xml"
fn = "data/BIOMD0000000469_url.xml" # big boi
d = sbml_to_sysinfo(fn)
td = make_nice(d)

rs = d["listOfRules"]
vars = map(x -> Symbol(x["variable"]), rs)

# todo find iv, make vars => vars(t)
dvs = @. Num(Variable(vars))
D = Differential(t)

eqs = @. parse_node(getproperty(rs, :firstelement))
eqs = reduce(vcat, eqs)

# iv problem, cant ODESystem(deqs)
deqs = D.(vars) .~ eqs 

pm = build_map(d["listOfParameters"], "id", "value")
sm = build_map(d["listOfSpecies"], "id", "initialConcentration")
cm = build_map(d["listOfCompartments"], "id", "size")


# do locally 
dfs = DataFrame.(values(td)) 
ddd = Dict(keys(td) .=> dfs)
