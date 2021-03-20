# This file uses MathML.jl and libsbml.jl to convert SBML models to ReactionSystems


mutable struct SbmlModel
    # Basic Functionality
    pars::Array{Pair{Sym,Float64}}  # [parameter=>value, ...]
    comps::Array{Pair{Sym,Float64}}  # [nucleus=>1.0,cytoplasm=>2.0,...]
    spec::Array{Num,Tuple{Float64,Num}}  # [specie=>(u0,compartment), ...]
    rxs::Array{Num,tuple{SymbolicUtils.Add,Array{Num},Array{Int},Array{Int}}} 
         # [(fullkineticlaw,[subtrate,...],[product,...],[substoich,...],[prodstoich,...]),...]
    ## Future components
    # maybe add fields to store info about events, piecewise simulation, etc.
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
    reactions = build_reactions(d["listOfReactions"])
    SbmlModel(pars,comps,spec,rxs)
end

function build_reactions(ps::Vector{EzXML.Node})
    reactions = []
    for reaction in ps
        reactants, r_stoich = _getlistofspeciesreference(reaction,"Reactant")
        products, p_stoich = _getlistofspeciesreference(reaction,"Product")
        kineticlaw = getkineticlaw(reaction)
        kineticlaw = parse_node(kineticlaw)
        thisreaction = (kineticlaw,reactants,products,r_stoich,p_stoich)
        append!(reactions,thisreaction)
    end
end

function _getlistofspeciesreference(reaction::EzXML.Node,type,name)
    if nodename(reaction) != "reaction"
        @error("Input node must be a reaction but is $(nodename(reaction)).")
    end
    listnodes = [node for node in eachelement(reaction) if nodename(node) == "listOf$(type)s"]
    if length(listnodes) > 1
        @error("SBML files with reactions with more than one listOf$(type)s are not supported.")
    end
    listnode = listnodes[1]
    spec = [Symbol(getindex(node, "id")) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    stoich = [getindex(node, "stoichiometry")) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    (spec, stoich)
end

function getkineticlaw(reaction::EzXML.Node)
    kineticlaws = [node for node in eachelement(reaction) if nodename(node) == "kineticLaw"]
    if length(kineticlaws) > 1
        @error("Reactions with more than one kinticLaw are not supported.")
    end
    kineticlaw = kineticlaws[1]
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
