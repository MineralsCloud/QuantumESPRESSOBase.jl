"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using MLStyle: @match
using Parameters: @with_kw
using Setfield: @lens, get

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards: Card,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    allowed_options,
    pseudopotential_format

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    KPoint,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    pseudopotential_format

# ================================== KPoint ================================== #
abstract type KPoint end

@with_kw struct MonkhorstPackGrid{A<:AbstractVector{Int},B<:AbstractVector{Int}} <: KPoint
    grid::A
    offsets::B
    @assert length(grid) == 3 "`grid` must be a three-element-vector! However it is of length $(length(grid))!"
    @assert length(offsets) == 3 "`offsets` must be a three-element-vector! However it is of length $(length(offsets))!"
    @assert all(x ∈ (0, 1) for x in offsets) "`offsets` must be either 0 or 1!"
end

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint{A<:AbstractVector{Float64},B<:Real} <: KPoint
    coordinates::A
    weight::B
    @assert length(coordinates) == 3 "`coordinates` must be a three-element-vector! However it is of length $(length(coordinates))!"
end
SpecialKPoint(x, y, z, w) = SpecialKPoint([x, y, z], w)

@with_kw struct KPointsCard{A<:AbstractString,B<:AbstractVector{<:KPoint}} <: Card
    option::A = "tpiba"
    data::B
    @assert option in allowed_options(KPointsCard)
    @assert begin
        @match option begin
            "automatic" => eltype(data) <: MonkhorstPackGrid
            "gamma" => eltype(data) <: GammaPoint
            # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            _ => eltype(data) <: SpecialKPoint
        end
    end
end
function KPointsCard(option::AbstractString, data::AbstractMatrix{<:Real})
    @assert(size(data, 2) == 4, "The size of the matrix should be (N, 4)! However $(size(data)) is given!")
    return KPointsCard(option, [SpecialKPoint(x...) for x in eachrow(data)])
end
# ============================================================================ #

allowed_options(::Type{<:KPointsCard}) =
    ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

end
