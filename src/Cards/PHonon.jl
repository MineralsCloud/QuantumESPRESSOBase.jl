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
using ..Cards

export QPoint, SpecialQPoint, QPointsSpecsCard

abstract type QPoint end

@with_kw struct SpecialQPoint{A<:AbstractVector{Float64},B<:Real} <: QPoint
    coordinates::A
    weight::B
    @assert length(coordinates) == 3
end

struct QPointsSpecsCard{A<:AbstractVector{SpecialQPoint}} <: Card
    data::A

@with_kw struct Q2RPoint{A<:AbstractVector{Float64}} <: QPoint
    coordinates::A
end

struct Q2RPointCard{A<:AbstractVector{Q2RPoint}} <: Card
    data::A
end

end
