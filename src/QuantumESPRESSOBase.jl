module QuantumESPRESSOBase

export InputEntry

abstract type InputEntry end

include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("bravais.jl")
include("Inputs/Inputs.jl")
include("name.jl")
include("to_qe.jl")
include("validate.jl")

end # module
