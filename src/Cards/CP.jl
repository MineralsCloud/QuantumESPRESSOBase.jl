module CP

using Rematch: @match
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
                                 potential_format

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
@with_kw struct AtomicVelocity{A<:AbstractVector{<:Real}}
    atom::String
    vel::A
    @assert(length(vel) == 3, "`vel` is not of length 3, but $(length(vel))!",)
end

struct AtomicVelocitiesCard{A<:AbstractVector{<:AtomicVelocity}} <: Card
    data::A
end
# ============================================================================ #

# ============================== RefCellParameters ============================== #
@with_kw struct RefCellParametersCard{A<:AbstractMatrix{<:Real}} <: AbstractCellParametersCard
    option::String = "bohr"
    data::A
    @assert(option ∈ allowed_options(RefCellParametersCard))
    @assert(size(data) == (3, 3))
end
# ============================================================================ #

# ============================== AtomicForce ============================== #
@with_kw struct AtomicForce{A<:AbstractVector{<:Real}}
    atom::String
    force::A
    @assert(length(force) == 3, "`force` is not of length 3, but $(length(force))!",)
end

@with_kw struct AtomicForcesCard{T<:AbstractVector{<:AtomicForce}} <: Card
    data::T
end
# ============================================================================ #

Cards.optionof(::AtomicVelocitiesCard) = "a.u"
Cards.optionof(::AtomicForcesCard) = nothing

Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

end
