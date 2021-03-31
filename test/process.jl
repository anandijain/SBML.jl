using SBML, Test, EzXML, ModelingToolkit
using CSV, DataFrames

cols = [:filename, :to_system, :to_problem, :to_solve, :states, :parameters, :error]

function test_sbml!(row, fn; dir="out")
	@show fn
    d, file  = splitdir(fn)
    b, ext = splitext(file)
    row[1] = fn
    try
        ml = ODESystem(fn)
        # sys = ml.sys
        row[2] = true
        row[5] = length(states(sys))
        row[6] = length(parameters(sys))
        prob = ODEProblem(ml, Pair[], (0., 1.))
        row[3] = true
        sol = solve(prob, Tsit5())
        row[4] = true
    catch e 
        row[end] = e
    end
    replace!(row, nothing => missing)
    CSV.write(joinpath(dir, file * ".txt"), DataFrame(cols .=> row))
    row
end

function test_sbmls!(mat, fns; dir="out")
    for i in eachindex(fns)
        @show i fns[i]
        mat[i, :] = test_sbml!(mat[i, :], fns[i];dir=dir)
    end
    mat
end

function test_sbml_suite(fns; dir="out")
    mkpath(dir)
    n = length(cols)
    mat = Array{Any,2}(nothing, length(fns), n)
    test_sbmls!(mat, fns; dir=dir)
    replace!(mat, nothing => missing)
    names .=> eachcol(mat) # used with DataFrame()
end

suite_dir = "../data/sbml"
biomodels_dir = "test/data/SBML_BIOMODELS"
@test isdir(suite_dir)
@test isdir(biomodels_dir)

suite_fns = readdir(suite_dir; join=true)
biomodels_fns = readdir(biomodels_dir; join=true)

test_sbml_suite(suite_fns; dir="out/suite/")
test_sbml_suite(biomodels_fns; dir="out/biomd/")

df = vcat(CSV.read.(readdir("out/suite/";join=true), DataFrame)...)
@test df isa DataFrame
CSV.write("sbml_test_suite_$(time()).csv", df)

dropmissing(df, :to_system)

df2 = vcat(CSV.read.(readdir("out/biomd/";join=true), DataFrame)...)
CSV.write("biomd_$(now()).csv", df)
