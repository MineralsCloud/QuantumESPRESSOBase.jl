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

export PhononInput

@with_kw struct PhononInput <: AbstractInput
    phonon::PhononNamelist = PhononNamelist()
    q_points::QPointsSpecsCard
    @assert phonon.qplot == true
end  # struct PhononInput

end
