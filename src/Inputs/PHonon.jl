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
using QuantumESPRESSOBase.Cards.PHonon: QPointsSpecsCard
using QuantumESPRESSOBase.Inputs: QuantumESPRESSOInput

export PhInput, Q2rInput, MatdynInput, DynmatInput

@with_kw struct PhInput <: QuantumESPRESSOInput
    inputph::PhNamelist = PhNamelist()
    q_points::Union{Nothing,QPointsSpecsCard} = nothing
end # struct PhInput

@with_kw struct Q2rInput <: QuantumESPRESSOInput
    input::Q2rNamelist = Q2rNamelist()
end

@with_kw struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist = MatdynNamelist()
    q_points::QPointsSpecsCard
end

@with_kw struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist = DynmatNamelist()
end

end
