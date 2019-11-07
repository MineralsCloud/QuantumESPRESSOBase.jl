"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase.Cards.PHonon
using QuantumESPRESSOBase.Inputs: QuantumESPRESSOInput

export PHononInput, Q2RInput, MatdynInput, DynmatInput

@with_kw struct PHononInput <: QuantumESPRESSOInput
    inputph::PHNamelist = PHNamelist()
    q_points::Union{Nothing,QPointsSpecsCard} = nothing
end # struct PHononInput

@with_kw struct Q2RInput <: QuantumESPRESSOInput
    input::Q2RNamelist = Q2RNamelist()
end

@with_kw struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist = MatdynNamelist()
    q_points::QPointsSpecsCard
end

@with_kw struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist = DynmatNamelist()
end

end
