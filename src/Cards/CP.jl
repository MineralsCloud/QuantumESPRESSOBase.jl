module CP

export UnifiedPseudopotentialFormat,
    VanderbiltUltraSoft,
    AndreaDalCorso,
    OldNormConserving,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    RefCellParametersCard,
    AtomicVelocity,
    AtomicVelocitiesCard,
    AtomicForce,
    AtomicForcesCard
export pseudopot_format, option_convert, push_atom!, append_atom!

include("shared.jl")

@auto_hash_equals mutable struct AtomicVelocity
    atom::String
    velocity::Vector{Float64}
    function AtomicVelocity(atom, velocity)
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        @assert length(velocity) == 3
        return new(atom, velocity)
    end
    function AtomicVelocity(atom::Union{AbstractChar,AbstractString})
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        return new(string(atom))
    end
end
AtomicVelocity(x::AtomicSpecies, velocity) = AtomicVelocity(x.atom, velocity)
AtomicVelocity(x::AtomicPosition, velocity) = AtomicVelocity(x.atom, velocity)
AtomicVelocity(x::AtomicSpecies) = AtomicVelocity(x.atom)
AtomicVelocity(x::AtomicPosition) = AtomicVelocity(x.atom)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicVelocity, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)
AtomicSpecies(x::AtomicVelocity) = AtomicSpecies(x.atom)
AtomicPosition(x::AtomicVelocity, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
AtomicPosition(x::AtomicVelocity) = AtomicPosition(x.atom)

@auto_hash_equals struct AtomicVelocitiesCard{A<:AbstractVector{<:AtomicVelocity}} <: Card
    data::A
end

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

Cards.optionof(::AtomicVelocitiesCard) = "a.u"
Cards.optionof(::AtomicForcesCard) = nothing

Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

function Base.setproperty!(value::AtomicVelocity, name::Symbol, x)
    if name == :atom
        @assert(length(x) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        if x isa AbstractChar
            x = string(x)
        end
    elseif name == :velocity && x isa AbstractVector  # It cannot be an `else` here, since it will capture invalid `name`s.
        @assert length(x) == 3
    end
    setfield!(value, name, x)
end # function Base.setproperty!

end
