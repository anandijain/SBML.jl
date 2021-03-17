using SBML
using Test

@testset "SBML.jl" begin
    # @testset "suite" begin include("test.jl") end
    @testset "data" begin include("data.jl") end
end
