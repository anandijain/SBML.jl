using MathML, EzXML, Symbolics

# const ns = namespace(node)
const sbml_lists = [
"listOfFunctionDefinitions",
"listOfUnitDefinitions",
"listOfCompartments",
"listOfSpecies",
"listOfParameters",
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
        str = "//x:$(list_name)"
        l = findall(str, node, ["x" => ns])
        #use children or elements?
        ls = !isempty(l) ? mapreduce(elements, vcat, l) : nothing
        push!(arr, list_name => ls)
    end
    Dict(arr)
end

