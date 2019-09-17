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
using QuantumESPRESSOBase.Inputs

export PHononInput, Q2RInput, MatdynInput, DynmatInput

@with_kw struct PHononInput <: AbstractInput
    inputph::PHNamelist = PHNamelist()
    q_points::QPointsSpecsCard
end # struct PHononInput

@with_kw struct Q2RInput <: AbstractInput
    input::Q2RNamelist = Q2RNamelist()
    grid::AbstractVector
    nqs::Int
    @assert length(grid) == 3
end

@with_kw struct MatdynInput <: AbstractInput
    input::MatdynNamelist = MatdynNamelist()
    q_points::QPointsSpecsCard
end

@with_kw struct DynmatInput <: AbstractInput
    input::DynmatNamelist = DynmatNamelist()
end

end
