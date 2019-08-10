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

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards

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

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies{A <: AbstractString,B <: Real,C <: AbstractString}
    atom::A
    mass::B
    pseudopotential::C
end

"""
    pseudopotential_format(data::AtomicSpecies)::String

Return the pseudopotential format.

The pseudopotential file is assumed to be in the new UPF format.
If it doesn't work, the pseudopotential format is determined by
the file name:
- "*.vdb or *.van": Vanderbilt US pseudopotential code
- "*.RRKJ3": Andrea Dal Corso's code (old format)
- none of the above: old PWscf norm-conserving format
"""
function pseudopotential_format(data::AtomicSpecies)::String
    @match lowercase(splitext(data.pseudopotential)[2]) begin
        "vdb" || "van" => "Vanderbilt US pseudopotential code"
        "rrkj3" => "Andrea Dal Corso's code (old format)"
        _ => "old PWscf norm-conserving format"
    end
end

struct AtomicSpeciesCard{T <: AbstractVector{<: AtomicSpecies}} <: Card
    data::T
end
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{A <: AbstractString,B <: AbstractVector{<: Real},C <: AbstractVector{Int}}
    atom::A
    pos::B; @assert length(pos) == 3
    if_pos::C = [1, 1, 1]; @assert length(if_pos) == 3
end

@with_kw struct AtomicPositionsCard{A <: AbstractString,B <: AbstractVector{<: AtomicPosition}} <: Card
    option::A = "alat"; @assert option in allowed_options(AtomicPositionsCard)
    data::B
end

function validate(x::AtomicSpeciesCard, y::AtomicPositionsCard)
    lens = @lens _.data.atom
    @assert isempty(symdiff(map(Base.Fix2(get, lens) ∘ unique, (x, y)))) "labels of the atoms are different in `ATOMIC_SPECIES` and `ATOMIC_POSITIONS` card!"
end # function validate
validate(y::AtomicPositionsCard, x::AtomicSpeciesCard) = validate(x, y)
# ============================================================================ #

# ============================== CellParameters ============================== #
@with_kw struct CellParametersCard{A <: AbstractString,B <: AbstractMatrix} <: Card
    option::A = "alat"; @assert option in allowed_options(CellParametersCard)
    data::B; @assert size(data) == (3, 3)
end
# ============================================================================ #

# ================================== KPoint ================================== #
abstract type KPoint end

@with_kw struct MonkhorstPackGrid{A <: AbstractVector{Int},B <: AbstractVector{Int}} <: KPoint
    grid::A; @assert length(grid) == 3
    offsets::B; @assert length(offsets) == 3 && all(x ∈ (0, 1) for x in offsets)
end

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint{A <: AbstractVector{Float64},B <: Real} <: KPoint
    coordinates::A; @assert length(coordinates) == 3
    weight::B
end

@with_kw struct KPointsCard{A <: AbstractString,B <: AbstractVector{<: KPoint}} <: Card
    option::A = "tpiba"; @assert option in allowed_options(KPointsCard)
    data::B
    @assert begin
        @match option begin
            "automatic" => eltype(data) <: MonkhorstPackGrid
            "gamma" => eltype(data) <: GammaPoint
            # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            _ => eltype(data) <: SpecialKPoint
        end
    end
end
# ============================================================================ #

# ================================== Methods ================================= #
eachrow(A::AbstractVecOrMat) = (view(A, i, :) for i in axes(A, 1))  # Julia 1.0 does not support `eachrow`
# ============================================================================ #

end
