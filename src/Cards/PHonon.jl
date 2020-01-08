"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Parameters: @with_kw

using QuantumESPRESSOBase.Cards: Card

import QuantumESPRESSOBase

export QPoint, SpecialQPoint, QPointsSpecsCard

abstract type QPoint end

@with_kw struct SpecialQPoint{A<:AbstractVector{Float64},B<:Real} <: QPoint
    coord::A
    weight::B
    @assert length(coord) == 3
end

struct QPointsSpecsCard{A<:AbstractVector{<:SpecialQPoint}} <: Card
    data::A
end

function QuantumESPRESSOBase.to_qe(
    card::QPointsSpecsCard;
    indent = ' '^4,
    delim = ' ',
)::String
    content = "$(length(card.data))\n"
    for p in card.data
        content *= indent * join([p.coord; p.weight], delim) * "\n"
    end
    return content
end

end
