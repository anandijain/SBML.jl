using SBML
using EzXML
using ModelingToolkit
using OrdinaryDiffEq
using Symbolics
using Test

@testset "SBML.jl" begin
    # @testset "suite" begin include("test.jl") end
    #=@testset "data" begin include("data.jl") end=#
    @testset "sbml2reactionsystem" begin include("sbml2reactionsystem.jl") end
end
