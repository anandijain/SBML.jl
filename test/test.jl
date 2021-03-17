using Test, SBML, EzXML, MathML, Symbolics

""" test on the sbml-test-suite, requires a clone of repo"""
function test_suite()
    folder = "C:\\Users\\Anand\\.julia\\dev\\data\\sbml"
    fns = readdir(folder; join=true)
    # 18 are bad
    arr = []
    for fn in fns 
        try
            d = sbml_to_sysinfo(fn)
            push!(arr, fn=>d)
        catch e end
    end
    @test length(fns) - length(arr) == 18
end

fn = "data/01193-sbml-l3v1.xml"
d = sbml_to_sysinfo(fn)
@test isa(d, Dict{String, Union{Nothing, Vector{EzXML.Node}}})
@test !all(isnothing, values(d))

fn = "data/BIOMD0000000346_url.xml"
d = sbml_to_sysinfo(fn)

rs = d["listOfRules"]
vars = map(x->Symbol(x["variable"]), rs)

params = d["listOfParameters"]
ss = d["listOfSpecies"]


"doesn't handle the `constant=Bool` attribute or units"
function build_parammap(ps::Vector{EzXML.Node})
    names = @. Symbol(getindex(ps, "name"))
    vals = @. Meta.parse(getindex(ps, "value"))
    names .=> vals
end
pm = build_parammap(params)

"doesn't handle the compartments"
function build_speciesmap(ps::Vector{EzXML.Node})
    names = @. Symbol(getindex(ps, "name"))
    u0 = @. Meta.parse(getindex(ps, "initialConcentration"))
    names .=> u0
end
sm = build_speciesmap(ss)


@. Num(Variable(vars))


@variables t x y
D = Differential(t)
eqs = @. parse_node(getproperty(rs, :firstelement))
eqvec = reduce(vcat, eqs)

deqs = D.(vars) .~ eqvec # lol this is probably wrong what is a "rate rule"
