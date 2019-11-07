module QuantumESPRESSOBase

abstract type InputEntry end

include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("bravais_lattice.jl")
include("Inputs/Inputs.jl")
include("postlude.jl")

end # module
