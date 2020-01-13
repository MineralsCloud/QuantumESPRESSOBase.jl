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
"""
    AtomicSpecies(x::AtomicVelocity, mass, pseudopot)
    AtomicSpecies(x::AtomicVelocity)

Construct an incomplete `AtomicSpecies` from an `AtomicVelocity` instance.
"""
AtomicSpecies(x::AtomicVelocity, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)
AtomicSpecies(x::AtomicVelocity) = AtomicSpecies(x.atom)
"""
    AtomicPosition(x::AtomicVelocity, pos, if_pos)
    AtomicPosition(x::AtomicVelocity)

Construct an incomplete `AtomicPosition` from an `AtomicVelocity` instance.
"""
AtomicPosition(x::AtomicVelocity, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
AtomicPosition(x::AtomicVelocity) = AtomicPosition(x.atom)

@auto_hash_equals struct AtomicVelocitiesCard{A<:AbstractVector{<:AtomicVelocity}} <: Card
    data::A
end

"""
    push_atom!(v::AbstractVector{AtomicVelocity}, atoms::AbstractString...)

Push an atom or multiple atoms to a vector of `AtomicVelocity`.

**Note**: these new `atoms` will result in incomplete `AtomicVelocity`s!

See also: [`push!`](@ref), [`append_atom!`](@ref)
"""
function push_atom!(v::AbstractVector{AtomicVelocity}, atoms::AbstractString...)
    return push!(v, map(AtomicVelocity, atoms)...)
end # function push_atom!
"""
    push_atom!(card::AtomicVelocitiesCard, atoms::AbstractString...)

Push an atom or multiple atoms to a `AtomicVelocitiesCard`.

**Note**: these new `atoms` will result in incomplete `AtomicVelocity`s!

See also: [`push!`](@ref), [`append_atom!`](@ref)
"""
function push_atom!(card::AtomicVelocitiesCard, atoms::AbstractString...)
    push!(card.data, map(AtomicVelocity, atoms)...)
    return card
end # function push_atom!

function append_atom!(
    v::AbstractVector{AtomicVelocity},
    atoms::AbstractVector{<:AbstractString},
)
    return append!(v, map(AtomicVelocity, atoms))
end # function append_atom!
function append_atom!(card::AtomicVelocitiesCard, atoms::AbstractVector{<:AbstractString})
    append!(card.data, map(AtomicVelocity, atoms))
    return card
end # function append_atom!

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
