# This file uses MathML.jl and libsbml.jl to convert SBML models to ReactionSystems


mutable struct SbmlModel
    # Basic Functionality
    parameters#::Array{Pair{Variable,Float64}}  # [parameter=>value, ...]
    compartments#::Array{Pair{Variable,Float64}}  # [nucleus=>1.0,cytoplasm=>2.0,...]
    species#::Array{Pair{Variable,Tuple{Float64,Variable}}}  # [specie=>(u0,compartment), ...]
    reactions#::Array{tuple{SymbolicUtils.Add,Array{Variable},Array{Variable},Array{Int},Array{Int}}} 
         # [(fullkineticlaw,[subtrate,...],[product,...],[substoich,...],[prodstoich,...]),...]
    ## Future components
    # maybe add fields to store info about events, piecewise simulation, conversionFactors etc.
end
SbmlModel(file::String) = SbmlModel(readxml(file))
SbmlModel(doc::EzXML.Document) = _process_doc(doc)


function _process_doc(doc)
    # doc = make_extensive(doc)
    doc = promotelocalparameters(doc)
    sysinfo = sbml_to_sysinfo(doc)
    pars = vcat(sysinfo["listOfParameters"],
                sysinfo["listOfReactions/x:reaction/x:kineticLaw/x:listOfParameters"],
                sysinfo["listOfReactions/x:reaction/x:kineticLaw/x:listOfLocalParameters"])
    pars = build_par_map(pars)

    comps = build_comp_map(sysinfo["listOfCompartments"])

    spec = build_spec_map(sysinfo["listOfSpecies"])
    reactions = build_reactions(sysinfo["listOfReactions"])
    SbmlModel(pars,comps,spec,reactions)
end

function build_reactions(listofreactions::Vector{EzXML.Node})
    reactions = Tuple{Num,Array{Num,1},Array{Num,1},Array{Int64,1},Array{Int64,1}}[]
    for reaction in listofreactions
        reactants, r_stoich = _getlistofspeciesreference(reaction,"Reactant")
        products, p_stoich = _getlistofspeciesreference(reaction,"Product")
        kineticlaw = getkineticlaw(reaction)
        kineticlaw = parse_node(getmath(kineticlaw))[1]
        println(typeof(reactants))
        println(typeof(r_stoich))
        println(typeof(products))
        println(typeof(p_stoich))
        println(typeof(kineticlaw))
        thisreaction = [Tuple{Num,Array{Num,1},Array{Num,1},Array{Int64,1},Array{Int64,1}}((kineticlaw,reactants,products,r_stoich,p_stoich))]
        println(typeof(reactions))
        println(typeof(thisreaction))
        append!(reactions,thisreaction)
    end
    reactions
end

function _getlistofspeciesreference(reaction::EzXML.Node,type)
    if nodename(reaction) != "reaction"
        @error("Input node must be a reaction but is $(nodename(reaction)).")
    end
    listnodes = [node for node in eachelement(reaction) if nodename(node) == "listOf$(type)s"]
    if length(listnodes) > 1
        @error("SBML files with reactions with more than one listOf$(type)s are not supported.")
    end
    listnode = listnodes[1]
    spec = [Num(Variable{Float64}(Symbol(getindex(node, "species")))) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    stoich = [Int(Meta.parse(getindex(node, "stoichiometry"))) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    (spec, stoich)
end

function getkineticlaw(reaction::EzXML.Node)
    kineticlaws = [node for node in eachelement(reaction) if nodename(node) == "kineticLaw"]
    if length(kineticlaws) > 1
        @error("Reactions with more than one kinticLaw are not supported.")
    end
    kineticlaw = kineticlaws[1]
end

function getmath(kineticlaw::EzXML.Node)
    maths = [node for node in eachelement(kineticlaw) if nodename(node) == "math"]
    if length(maths) > 1
        @error("kineticLaws with more than one math node are not supported.")
    end
    math = maths[1]
end

function build_par_map(parnodes::Vector{EzXML.Node})
    ids = @. Num(Variable{Float64}(Symbol(getindex(parnodes, "id"))))
    vals = @. Float64(Meta.parse(getindex(parnodes, "value")))
    ids .=> vals
end

function build_comp_map(compnodes::Vector{EzXML.Node})
    ids = @. Num(Variable{Float64}(Symbol(getindex(compnodes, "id"))))
    sizes = @. Float64(Meta.parse(getindex(compnodes, "size")))
    ids .=> sizes
end

function build_spec_map(compnodes::Vector{EzXML.Node})
    ids = @. Num(Variable{Float64}(Symbol(getindex(compnodes, "id"))))
    inits = @. Float64(Meta.parse(getindex(compnodes, "initialAmount")))
    comps = @. Num(Variable{Float64}(Symbol(getindex(compnodes, "compartment"))))
    ids .=> tuple.(inits,comps)
end

#=function expand_reversible(doc::EzXML.Node)
    doc = deepcopy(doc)
    root = doc.root
    ns = namespace(root)
    reactions = findall("//x:reaction", root, ["x"=>ns])
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
            par["name"] = reaction["name"]*"_"*par["name"]
        end
    end
    doc
end=#

function promotelocalparameters(doc::EzXML.Document)
    doc = deepcopy(doc)
    root = doc.root
    ns = namespace(root)
    reactions = findall("//x:reaction", root, ["x"=>ns])
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
            par["name"] = reaction["name"]*"_"*par["name"]
        end
    end
    doc
end
