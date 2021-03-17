using Test, SBML, EzXML, MathML, Symbolics

fn = "data/01193-sbml-l3v1.xml"
d = sbml_to_sysinfo(fn)
@test isa(d, Dict{String, Union{Nothing, Vector{EzXML.Node}}})
@test !all(isnothing, values(d))

fn = "data/BIOMD0000000346_url.xml"
d = sbml_to_sysinfo(fn)
ks, vs = keys(d), values(d)
ts = to_table.(values(d))
td = Dict(ks .=> ts)
filter!(x->!isnothing(x.second), td)

@test isa(td, Dict)

# Namespacing issue: fbc:ChemicalFormula
@test isempty(collect(skipmissing(td["listOfSpecies"][end].second)))

# do locally, i don't want DF dep
# dfs = @. DataFrame(last(td)) 
# ddd = Dict(keys(td) .=> dfs)
