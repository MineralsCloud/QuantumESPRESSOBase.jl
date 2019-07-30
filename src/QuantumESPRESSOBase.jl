module QuantumESPRESSOBase

export InputEntry

abstract type InputEntry end

include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("QuantumESPRESSOInput/QuantumESPRESSOInput.jl")

end # module
