module CP

using MLStyle: @match
using Parameters: @with_kw

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards: Card,
                                 AtomicSpecies,
                                 AtomicSpeciesCard,
                                 AtomicPosition,
                                 AtomicPositionsCard,
                                 AbstractCellParametersCard,
                                 CellParametersCard,
                                 pseudopotential_format

export AtomicSpecies,
       AtomicSpeciesCard,
       AtomicPosition,
       AtomicPositionsCard,
       CellParametersCard,
       AtomicVelocity,
       AtomicVelocitiesCard,
       RefCellParametersCard,
       AtomicForce,
       AtomicForcesCard

# ============================== AtomicVelocity ============================== #
@with_kw struct AtomicVelocity{A<:AbstractString,B<:AbstractVector{<:Real}}
    atom::A
    vel::B
    @assert(
        length(vel) == 3,
        "`vel` must be a three-element-vector! However it is of length $(length(vel))!",
    )
end

struct AtomicVelocitiesCard{A<:AbstractVector{<:AtomicVelocity}} <: Card
    data::A
end
# ============================================================================ #

# ============================== RefCellParameters ============================== #
@with_kw struct RefCellParametersCard{
    A<:AbstractString,
    B<:AbstractMatrix,
} <: AbstractCellParametersCard
    option::A = "bohr"
    data::B
    @assert(option ∈ allowed_options(RefCellParametersCard))
    @assert(size(data) == (3, 3))
end
# ============================================================================ #

# ============================== AtomicForce ============================== #
@with_kw struct AtomicForce{A<:AbstractString,B<:AbstractVector{<:Real}}
    atom::A
    force::B
    @assert(
        length(force) == 3,
        "`force` must be a three-element-vector! However it is of length $(length(force))!",
    )
end

@with_kw struct AtomicForcesCard{T<:AbstractVector{<:AtomicForce}} <: Card
    data::T
end
# ============================================================================ #

Cards.option(::AtomicVelocitiesCard) = "a.u"
Cards.option(::AtomicForcesCard) = nothing

Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

const ANGSTROM_TO_BOHR = 1 / 0.529177210903

"""
    cell_volume(card::RefCellParametersCard)

Return the cell volume according to the `RefCellParametersCard`'s parameters, in atomic unit.
"""
function Cards.cell_volume(card::RefCellParametersCard)
    @match option(card) begin
        "bohr" => det(card.data)
        "angstrom" => det(card.data) * ANGSTROM_TO_BOHR^3
    end
end # function cell_volume

end
