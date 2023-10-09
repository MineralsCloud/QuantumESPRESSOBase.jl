module PHonon

using Setfield: @set!
using StructEquality: @struct_hash_equal

using ..QuantumESPRESSOBase: ReciprocalPoint, Card

import ..QuantumESPRESSOBase: VerbositySetter, groupname

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

@struct_hash_equal struct QPointsCard <: Card
    data::Vector{ReciprocalPoint}
end
function QPointsCard(data::AbstractMatrix)
    @assert size(data, 2) == 4
    return QPointsCard(map(ReciprocalPoint, eachrow(data)))
end

include("inputs.jl")

end
