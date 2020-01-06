module CP

using Parameters: @with_kw

using QuantumESPRESSOBase.Cards: Card,
                                 AtomicSpecies,
                                 AtomicSpeciesCard,
                                 AtomicPosition,
                                 AtomicPositionsCard,
                                 AbstractCellParametersCard,
                                 CellParametersCard,
                                 AtomicForce,
                                 AtomicForcesCard

import QuantumESPRESSOBase.Cards

export AtomicSpecies,
       AtomicSpeciesCard,
       AtomicPosition,
       AtomicPositionsCard,
       CellParametersCard,
       AtomicForce,
       AtomicForcesCard,
       AtomicVelocity,
       AtomicVelocitiesCard,
       RefCellParametersCard,
       AtomicForce,
       AtomicForcesCard

# ============================== AtomicVelocity ============================== #
struct AtomicVelocity{A<:AbstractVector{<:Real}}
    atom::String
    v::A
    function AtomicVelocity{A}(atom, v) where {A<:AbstractVector{<:Real}}
        @assert length(v) == 3
        return new(atom, v)
    end # function AtomicVelocity
end
AtomicVelocity(atom, v::A) where {A} = AtomicVelocity{A}(atom, v)

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

Cards.optionof(::AtomicVelocitiesCard) = "a.u"
Cards.optionof(::AtomicForcesCard) = nothing

Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

end
