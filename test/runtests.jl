using QuantumESPRESSOBase
using Test

@testset "QuantumESPRESSOBase.jl" begin
    include("Inputs/PWscf.jl")
    include("Cards/Cards.jl")
    include("Setters.jl")
end
