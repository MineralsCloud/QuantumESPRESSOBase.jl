using QuantumESPRESSOBase
using Test

@testset "QuantumESPRESSOBase.jl" begin
    include("PWscf.jl")
    include("PHonon.jl")
    include("set.jl")
    # include("CP.jl")
end
