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

export PHononInput, Q2RInput, MatdynInput

@with_kw struct PHononInput <: AbstractInput
    phonon::PhononNamelist = PhononNamelist()
    q_points::QPointsSpecsCard
    @assert phonon.qplot == true
end # struct PHononInput

@with_kw struct Q2RInput <: AbstractInput
    q2r::Q2RNamelist = Q2RNamelist()
    grid::AbstractVector
    nqs::Int
    @assert length(grid) == 3
end

@with_kw struct MatdynInput <: AbstractInput
    matdyn::MatdynNamelist = MatdynNamelist()
    q_points::QPointsSpecsCard
end

end
