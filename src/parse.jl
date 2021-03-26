using MathML, EzXML, Symbolics

# const ns = namespace(node)
const sbml_lists = [
"listOfFunctionDefinitions",
"listOfUnitDefinitions",
"listOfCompartments",
"listOfSpecies",
"listOfParameters",
"listOfReactions/x:reaction/x:listOfParameters",
"listOfReactions/x:reaction/x:listOfLocalParameters",
"listOfInitialAssignments",
"listOfRules",
"listOfConstraints",
"listOfReactions",
"listOfEvents",
]

"""
give a fn, doc, or node of an sbml model and gather metadata

this gets used to build the MTK component
"""
function sbml_to_sysinfo end

function sbml_to_sysinfo(fn::AbstractString)
    !isfile(fn) && error("isfile(fn) must be true")
    sbml_to_sysinfo(readxml(fn))
end

function sbml_to_sysinfo(doc::EzXML.Document)
    sbml_to_sysinfo(doc.root)
end

function sbml_to_sysinfo(node::EzXML.Node)
    ns = namespace(node) # will this always work
    node.name != "sbml" && error("give an sbml node")
    arr = []
    for list_name in sbml_lists 
        str = "/x:sbml/x:model/x:$(list_name)"
        l = findall(str, node, ["x" => ns])
        #use children or elements?
        ls = !isempty(l) ? mapreduce(elements, vcat, l) : nothing
        push!(arr, list_name => ls)
    end
    Dict(arr)
end

#=function sbml_to_sysinfo(node::EzXML.Node,promotelocalparameters::Bool)
    ns = namespace(node) # will this always work
    node.name != "sbml" && error("give an sbml node")
    arr = []
    for list_name in sbml_lists 
        str = "/x:$(list_name)"
        l = findall(str, node, ["x" => ns])
        #use children or elements?
        ls = !isempty(l) ? mapreduce(elements, vcat, l) : nothing
        push!(arr, list_name => ls)
    end
    Dict(arr)
end=#

"doesn't handle the `constant=Bool` attribute or units"
function build_map(ps::Vector{EzXML.Node}, name, value)
    names = @. Symbol(getindex(ps, name))
    vals = @. Meta.parse(getindex(ps, value))
    names .=> vals
end

"given a list of ezxml nodes, return the union of all unique attributes of "
function uniquekeys(node_list)
    isnothing(node_list) ? nothing :
    unique(get_name.(reduce(vcat, attributes.(node_list))))
end

"create a `DataFrame`-able for the node list, returns pairs"
function to_table(node_list)
    isnothing(node_list) && return nothing
    ks = uniquekeys(node_list)
    arr = []
    for node in node_list
        r = []
        for k in ks
            x = haskey(node, k) ? node[k] : missing
            push!(r, x)
        end
        push!(arr, r)
    end
    ks .=> eachrow(reduce(hcat, arr))
end

get_name(n::EzXML.Node) = n.name

# function build_tables(fn)
#     ts = Dict()
#     sysinfo = sbml_to_sysinfo(fn)
#     for (k, v) in sysinfo
#         @show k v
#         t = isnothing(v) ? to_table(v) : nothing   
#         ts[k] = t
#     end
#     ts
# end

function make_nice(d)
    ks, vs = keys(d), values(d)
    ts = to_table.(values(d))
    td = Dict(ks .=> ts)
    filter(x -> !isnothing(x.second), td)
end