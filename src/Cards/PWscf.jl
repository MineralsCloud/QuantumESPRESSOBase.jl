"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: eachrow
using Rematch: @match
using Parameters: @with_kw
using Setfield: @lens, get

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards: Card,
                                 AtomicSpecies,
                                 AtomicSpeciesCard,
                                 AtomicPosition,
                                 AtomicPositionsCard,
                                 CellParametersCard,
                                 AtomicForce,
                                 AtomicForcesCard,
                                 potential_format

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
       potential_format

# ================================== KPoint ================================== #
abstract type KPoint end

@with_kw struct MonkhorstPackGrid{A<:AbstractVector{<:Integer},B<:AbstractVector{<:Integer}}
    grid::A
    offsets::B
    @assert(length(grid) == 3, "`grid` is not of length 3, but $(length(grid))!",)
    @assert(length(offsets) == 3, "`offsets` is not of length 3, but $(length(offsets))!",)
    @assert(all(x ∈ (0, 1) for x in offsets), "`offsets` must be either 0 or 1!")
end

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint{A<:AbstractVector{<:Real},B<:Real} <: KPoint
    coordinates::A
    weight::B
    @assert(
        length(coordinates) == 3,
        "`coordinates` is not of length 3, but $(length(coordinates))!",
    )
end
SpecialKPoint(x, y, z, w) = SpecialKPoint([x, y, z], w)

@with_kw struct KPointsCard{A<:Union{
    MonkhorstPackGrid,
    GammaPoint,
    AbstractVector{<:SpecialKPoint},
}} <: Card
    option::String = "tpiba"
    data::A
    @assert(option ∈ allowed_options(KPointsCard))
    @assert begin
        @match option begin
            "automatic" => typeof(data) <: MonkhorstPackGrid
            "gamma" => typeof(data) <: GammaPoint
            # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            _ => eltype(data) <: SpecialKPoint
        end
    end
end
function KPointsCard(option::AbstractString, data::AbstractMatrix{<:Real})
    @assert(size(data, 2) == 4, "The size of `data` is not `(N, 4)`, but $(size(data))!",)
    return KPointsCard(option, [SpecialKPoint(x...) for x in eachrow(data)])
end
# ============================================================================ #

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

end
