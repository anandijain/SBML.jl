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

create_var(x) = Num(Variable(Symbol(x)))
# create_var(x, iv) = Num(Sym{FnType{Tuple{Real}}}(Symbol(x))(Variable(Symbol(iv)))).val
# create_var(x, iv) = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(Symbol(x)))(Variable(Symbol(iv)))
create_param(x) = Num(Sym{ModelingToolkit.Parameter{Real}}(Symbol(x)))


""" Convert SBML document into SbmlModel """
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
    substitutions = _mathml_substitutions(sysinfo)
    reactions = build_reactions(sysinfo["listOfReactions"], substitutions)
    SbmlModel(pars,comps,spec,reactions)
end

function _mathml_substitutions(sysinfo)
    pars = vcat(sysinfo["listOfCompartments"],
                sysinfo["listOfParameters"],
                sysinfo["listOfReactions/x:reaction/x:kineticLaw/x:listOfParameters"],
                sysinfo["listOfReactions/x:reaction/x:kineticLaw/x:listOfLocalParameters"])
    vars = sysinfo["listOfSpecies"]
    par_subs = _par_subs(pars)
    # var_subs = _var_subs(vars)
    Dict(par_subs)
end


function _par_subs(nodes)
    mathml = @. create_var(getindex(nodes, "id"))
    sbmlmodel = @. create_param(getindex(nodes, "id"))
    mathml .=> sbmlmodel
end


#=function _var_subs(nodes)
    mathml = @. Num(Variable(Symbol(getindex(nodes, "id"))))
    sbmlmodel = @. Num(Sym{ModelingToolkit.Parameter{Real}}(Symbol(getindex(nodes, "id"))))
    mathml .=> sbmlmodel
end=#


""" Extract kineticLaws, Reactants, Products and their stoichiometry from SBML """
function build_reactions(listofreactions::Vector{EzXML.Node}, substitutions)
    reactions = Tuple{Num,Array{Num,1},Array{Num,1},Array{Int64,1},Array{Int64,1}}[]
    for reaction in listofreactions
        reactants, r_stoich = _getlistofspeciesreference(reaction,"Reactant")
        products, p_stoich = _getlistofspeciesreference(reaction,"Product")
        kineticlaw = getkineticlaw(reaction)
        kineticlaw = parse_node(getmath(kineticlaw))[1]
        kineticlaw = substitute(kineticlaw, substitutions)
        thisreaction = [Tuple{Num,Array{Num,1},Array{Num,1},Array{Int64,1},Array{Int64,1}}((kineticlaw,reactants,products,r_stoich,p_stoich))]
        append!(reactions,thisreaction)
    end
    reactions
end

""" Extract species and their stroichiometry from an SBML reaction """
function _getlistofspeciesreference(reaction::EzXML.Node,type)
    if nodename(reaction) != "reaction"
        @error("Input node must be a reaction but is $(nodename(reaction)).")
    end
    listnodes = [node for node in eachelement(reaction) if nodename(node) == "listOf$(type)s"]
    if length(listnodes) > 1
        @error("SBML files with reactions with more than one listOf$(type)s are not supported.")
    end
    listnode = listnodes[1]
    spec = [create_var(getindex(node, "species")) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    stoich = [Int(Meta.parse(getindex(node, "stoichiometry"))) for node in eachelement(listnode) if nodename(node) == "speciesReference"]
    (spec, stoich)
end

""" Extract kineticLaw node from an SBML reaction """
function getkineticlaw(reaction::EzXML.Node)
    kineticlaws = [node for node in eachelement(reaction) if nodename(node) == "kineticLaw"]
    if length(kineticlaws) > 1
        @error("Reactions with more than one kinticLaw are not supported.")
    end
    kineticlaw = kineticlaws[1]
end

""" Extract math node from an SBML kineticLaw """
function getmath(kineticlaw::EzXML.Node)
    maths = [node for node in eachelement(kineticlaw) if nodename(node) == "math"]
    if length(maths) > 1
        @error("kineticLaws with more than one math node are not supported.")
    end
    math = maths[1]
end

""" Extract parameters from SBML """
function build_par_map(parnodes::Vector{EzXML.Node})
    ids = @. create_param(getindex(parnodes, "id"))
    vals = @. Float64(Meta.parse(getindex(parnodes, "value")))
    ids .=> vals
end

""" Extract compartments from SBML """
function build_comp_map(compnodes::Vector{EzXML.Node})
    ids = @. create_param(getindex(compnodes, "id"))
    sizes = @. Float64(Meta.parse(getindex(compnodes, "size")))
    ids .=> sizes
end

""" Extract species from SBML """
function build_spec_map(compnodes::Vector{EzXML.Node})
    # ids = @. Num(Variable{Float64}(Symbol(getindex(compnodes, "id"))))
    ids = @. create_var(getindex(compnodes, "id"))
    inits = @. Float64(Meta.parse(getindex(compnodes, "initialAmount")))
    comps = @. create_param(getindex(compnodes, "compartment"))
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

""" Assign gloablly unique ids to local parameters """
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
            _substitute_par!(kineticlaw, par["id"], reaction["id"])
            par["id"] = reaction["id"]*"_"*par["id"]
            par["name"] = reaction["id"]*"_"*par["name"]
        end
    end
    doc
end

function _substitute_par!(node, par, prefix)
    if typeof(node) != Nothing
        for element in eachelement(node)
            _substitute_par!(element, par, prefix)
        end
        if strip(node.name) == "ci" && strip(node.content) == par
            node.content = prefix*"_"*par
        end
    end
end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(sbmlmodel::SbmlModel)
    rxs = [Reaction(reac...; use_only_rate=true) for reac in sbmlmodel.reactions]
    t = Variable(:t)
    species = [spec.first for spec in sbmlmodel.species]
    pc = append!(sbmlmodel.parameters,sbmlmodel.compartments)
    params = [par.first for par in pc]
    ReactionSystem(rxs,t,species,params)
end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(sbmldocument::EzXML.Document)
    sbmlmodel = SbmlModel(sbmldocument)
    ReactionSystem(sbmlmodel)
end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(sbmlfile::String)
    sbmlmodel = SbmlModel(sbmlfile)
    ReactionSystem(sbmlmodel)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(sbmlmodel::SbmlModel)
    rs = ReactionSystem(sbmlmodel)
    odesys = convert(ODESystem, rs)
end

""" ODEProblem constructor """
function ModelingToolkit.ODEProblem(sbmlmodel::SbmlModel,tspan)
    odesys = ODESystem(sbmlmodel)
    u0map = [s.first => s.second[1] for s in sbmlmodel.species]
    parammap = append!(sbmlmodel.parameters,sbmlmodel.compartments)
    ODEProblem(odesys, u0map, tspan, parammap)
end
