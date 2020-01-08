"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using AutoHashEquals: @auto_hash_equals
using Compat: eachrow
using Formatting: sprintf1
using Parameters: @with_kw
using Setfield: @lens, get

using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Cards:
    Card,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    UnifiedPseudopotentialFormat,
    VanderbiltUltraSoft,
    AndreaDalCorso,
    OldNormConserving,
    pseudopot_format,
    allowed_options

import QuantumESPRESSOBase
import QuantumESPRESSOBase.Cards

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    KPoint,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    UnifiedPseudopotentialFormat,
    VanderbiltUltraSoft,
    AndreaDalCorso,
    OldNormConserving,
    pseudopot_format

abstract type KPoint end

@auto_hash_equals struct MonkhorstPackGrid{
    A<:AbstractVector{<:Integer},
    B<:AbstractVector{<:Integer},
}
    grid::A
    offsets::B
    function MonkhorstPackGrid{A,B}(grid, offsets) where {A,B}
        @assert(length(grid) == length(offsets) == 3)
        # See https://github.com/aiidateam/aiida-quantumespresso/blob/4aef9f9/aiida_quantumespresso/cli/utils/validate.py#L10-L37
        @assert(all(grid .> 0), "`grid` must be positive integers!")
        @assert(all(iszero(x) || isone(x) for x in offsets), "`offsets` must be 0 or 1!")
        return new(grid, offsets)
    end # function MonkhorstPackGrid
end
MonkhorstPackGrid(grid::A, offsets::B) where {A,B} = MonkhorstPackGrid{A,B}(grid, offsets)

struct GammaPoint <: KPoint end

@auto_hash_equals struct SpecialKPoint{A<:AbstractVector{<:Real},B<:Real} <: KPoint
    coord::A
    weight::B
    function SpecialKPoint{A,B}(coord, weight) where {A,B}
        @assert(length(coord) == 3)
        return new(coord, weight)
    end # function SpecialKPoint
end
SpecialKPoint(coord::A, weight::B) where {A,B} = SpecialKPoint{A,B}(coord, weight)
SpecialKPoint(x, y, z, w) = SpecialKPoint([x, y, z], w)

@with_kw struct KPointsCard{
    A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{<:SpecialKPoint}},
} <: Card
    option::String = "tpiba"
    data::A
    @assert(option ∈ allowed_options(KPointsCard))
    @assert if option == "automatic"
        typeof(data) <: MonkhorstPackGrid
    elseif option == "gamma"
        typeof(data) <: GammaPoint
    else  # option ∈ ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        eltype(data) <: SpecialKPoint
    end
end
function KPointsCard(option::AbstractString, data::AbstractMatrix{<:Real})
    @assert(size(data, 2) == 4, "The size of `data` is not `(N, 4)`, but $(size(data))!")
    return KPointsCard(option, [SpecialKPoint(x...) for x in eachrow(data)])
end

Cards.allowed_options(::Type{<:KPointsCard}) = (
    "tpiba",
    "automatic",
    "crystal",
    "gamma",
    "tpiba_b",
    "crystal_b",
    "tpiba_c",
    "crystal_c",
)

QuantumESPRESSOBase.asfieldname(::Type{<:KPointsCard}) = :k_points

QuantumESPRESSOBase.titleof(::Type{<:KPointsCard}) = "K_POINTS"

QuantumESPRESSOBase.to_qe(data::GammaPoint) = ""
function QuantumESPRESSOBase.to_qe(
    data::Union{MonkhorstPackGrid,SpecialKPoint};
    delim = ' ',
    numfmt = "%14.9f",
)
    return join(
        map(x -> sprintf1(numfmt, x), [getfield(data, 1); getfield(data, 2)]),
        delim,
    )
end
function QuantumESPRESSOBase.to_qe(
    card::KPointsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
)
    content = "K_POINTS { $(card.option) }\n"
    if optionof(card) in ("gamma", "automatic")
        content *= indent * to_qe(card.data) * "\n"
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= "$(length(card.data))\n"
        for x in card.data
            content *= indent * to_qe(x; delim = delim, numfmt = numfmt) * "\n"
        end
    end
    return content
end

end
