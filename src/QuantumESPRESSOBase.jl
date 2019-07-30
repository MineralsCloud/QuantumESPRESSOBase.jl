module QuantumESPRESSOBase

export InputEntry

abstract type InputEntry end

include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("QuantumESPRESSOInput/QuantumESPRESSOInput.jl")
include("name.jl")
include("to_qe.jl")

end # module
