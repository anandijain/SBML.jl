# This file uses MathML.jl and libsbml.jl to convert SBML models to ReactionSystems


mutable struct SbmlModel
    # Basic Functionality
    pars::Dict{Num,Float64}  # Dict(parameter=>value, ...)
    comps::Array{Num}  # [nucleus,cytoplasm,...]
    spec::Dict{Num,Tuple{Float64,Num}}  # Dict(specie=>(u0,compartment), ...)
    rxs::Array{Num,tuple{SymbolicUtils.Add,Array{Num},Array{Num},Array{Int},Array{Int}}} 
         # [(massactionrateconstant,fullkineticlaw,[subtrate,...],[product,...],[substoich,...],[prodstoich,...]),...]
    ## Future components
    # maybe add fields to store info about events, piecewise simulation, etc.
    function SbmlModel()
        return new(
            Dict(),
            Array
            Dict{String,T}(),
            []
        )
    end
end
SbmlModel(doc::XMLDocument) = process_doc(doc)

function process_doc(doc)
    doc = make_extensive(doc)
    doc = promotelocalparameters(doc)
    d = sbml_to_sysinfo(doc)
    globalpars = build_map(d["listOfParameters"], "id", "value")
    localpars = build_map(d["listOfLocalParameters"], "id", "value")
    pars = append!(globalpars, localpars)

    cm = build_map(d["listOfCompartments"], "id", "size")

    spec = build_map(d["listOfSpecies"], "id", "initialAmount", "compartment")

    # Todo: Reaction

end



function build_map(ps::Vector{EzXML.Node}, name, values...)
    names = @. Symbol(getindex(ps, name))
    vals = (@. Meta.parse(getindex(ps, value)) for value in values...)
    names .=> vals
end


function promotelocalparameters(doc::EzXML.Node)
    doc = deepcopy(doc)
    ns = namespace(doc)
    reactions = findall("//x:reaction", model, ["x"=>ns])
    for reaction in reactions
        kineticlaws = [node for node in eachelement(reaction) if nodename(node) == "kineticLaw"]
        if length(kineticlaws) > 1
            @error("SBML reactions with more than one kinetic law are not supported.")
        end
        kineticlaw = kineticlaws[1]
        locallistsofparameters = [node for node in eachelement(kineticlaw) if occursin("Parameters",nodename(node))]
        locallistofparameters = [node for list in locallistsofparameters for node in eachelement(list)]            
        for par in locallistofparameters
            par["id"] = reaction["id"]*"_"*par["id"]
        end
    doc
end
