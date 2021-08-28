"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using AbInitioSoftwareBase.Inputs: Namelist
using AutoHashEquals: @auto_hash_equals
using Compat: @NamedTuple
using ConstructionBase: setproperties
using Crystallography: ReciprocalPoint
using Setfield: @set!

using ..Inputs: Card, QuantumESPRESSOInput, VerbositySetter
using ..Inputs.PWscf: PWInput

import ..Inputs: groupname

export QPointsCard,
    PhInput,
    Q2rInput,
    MatdynInput,
    DynmatInput,
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist,
    VerbositySetter
export relayinfo

include("ph.jl")
include("q2r.jl")
include("matdyn.jl")
include("dynmat.jl")

groupname(::Type{PhNamelist}) = "INPUTPH"
groupname(::Type{Q2rNamelist}) = "INPUT"
groupname(::Type{MatdynNamelist}) = "INPUT"
groupname(::Type{DynmatNamelist}) = "INPUT"

struct QPointsCard <: Card
    data::Vector{ReciprocalPoint}
end

include("inputs.jl")

end
