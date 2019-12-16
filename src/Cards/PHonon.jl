"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Parameters: @with_kw

using QuantumESPRESSOBase
using ..Cards: Card

export QPoint, SpecialQPoint, QPointsSpecsCard

abstract type QPoint end

@with_kw struct SpecialQPoint{A<:AbstractVector{Float64},B<:Real} <: QPoint
    coordinates::A
    weight::B
    @assert length(coordinates) == 3
end

struct QPointsSpecsCard{A<:AbstractVector{<:SpecialQPoint}} <: Card
    data::A
end

function QuantumESPRESSOBase.to_qe(
    card::QPointsSpecsCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)::String
    content = "$(length(card.data))\n"
    for p in card.data
        content *= indent * join([p.coordinates; p.weight], sep) * "\n"
    end
    return content
end

end
