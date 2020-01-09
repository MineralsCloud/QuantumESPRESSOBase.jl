"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist
using QuantumESPRESSOBase.Cards.PHonon: KPointsCard
using QuantumESPRESSOBase.Inputs: QuantumESPRESSOInput

export PhInput, Q2rInput, MatdynInput, DynmatInput

@with_kw struct PhInput <: QuantumESPRESSOInput
    inputph::PhNamelist = PhNamelist()
    q_points::Union{Nothing,KPointsCard} = nothing
end # struct PhInput

struct Q2rInput <: QuantumESPRESSOInput
    input::Q2rNamelist
end # struct Q2rInput
Q2rInput() = Q2rInput(Q2rNamelist())

@with_kw struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist = MatdynNamelist()
    q_points::KPointsCard
end # struct MatdynInput

struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist
end # struct DynmatInput
DynmatInput() = DynmatInput(DynmatNamelist())

end
