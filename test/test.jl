using Test, SBML, EzXML, MathML, Symbolics

# folder = "C:\\Users\\Anand\\.julia\\dev\\data\\sbml"

""" test on the sbml-test-suite, requires a clone of repo"""
function test_suite(cases)
    fns = readdir(cases; join=true)
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
