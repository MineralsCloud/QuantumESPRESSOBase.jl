"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Crystallography: ReciprocalPoint
using Setfield: @set!

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

include("inputs.jl")

end
