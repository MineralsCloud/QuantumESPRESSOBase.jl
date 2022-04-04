"""
# module PHonon



# Examples

```julia
julia>
```
"""
module PHonon

using Crystallography: ReciprocalPoint
using Setfield: @set!
using StructHelpers: @batteries

using ..Inputs: Card

import ..Inputs: VerbositySetter, groupname

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

include("namelists.jl")

struct QPointsCard <: Card
    data::Vector{ReciprocalPoint}
end
@batteries QPointsCard eq = true hash = true

include("inputs.jl")

end
