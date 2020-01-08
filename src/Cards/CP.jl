module CP

using AutoHashEquals: @auto_hash_equals

using QuantumESPRESSOBase.Cards:
    Card,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    AbstractCellParametersCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    UnifiedPseudopotentialFormat,
    VanderbiltUltraSoft,
    AndreaDalCorso,
    OldNormConserving,
    pseudopot_format

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
    AtomicForcesCard,
    UnifiedPseudopotentialFormat,
    VanderbiltUltraSoft,
    AndreaDalCorso,
    OldNormConserving,
    pseudopot_format

# ============================== AtomicVelocity ============================== #
@auto_hash_equals struct AtomicVelocity{A<:AbstractVector{<:Real}}
    atom::String
    velocity::A
    function AtomicVelocity{A}(atom, velocity) where {A<:AbstractVector{<:Real}}
        @assert length(velocity) == 3
        return new(atom, velocity)
    end
end
AtomicVelocity(atom, velocity::A) where {A} = AtomicVelocity{A}(atom, velocity)

@auto_hash_equals struct AtomicVelocitiesCard{A<:AbstractVector{<:AtomicVelocity}} <: Card
    data::A
end
# ============================================================================ #

# ============================== RefCellParameters ============================== #
@auto_hash_equals struct RefCellParametersCard{A<:AbstractMatrix{<:Real}} <:
                         AbstractCellParametersCard
    option::String
    data::A
    function RefCellParametersCard{A}(option, data) where {A<:AbstractMatrix{<:Real}}
        @assert option ∈ allowed_options(RefCellParametersCard)
        @assert size(data) == (3, 3)
        return new(option, data)
    end
end
RefCellParametersCard(option, data::A) where {A} = RefCellParametersCard{A}(option, data)
RefCellParametersCard(data) = RefCellParametersCard("bohr", data)
# ============================================================================ #

Cards.optionof(::AtomicVelocitiesCard) = "a.u"
Cards.optionof(::AtomicForcesCard) = nothing

Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

end
