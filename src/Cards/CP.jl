module CP

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    RefCellParametersCard,
    AtomicVelocity,
    AtomicVelocitiesCard,
    AtomicForce,
    AtomicForcesCard
export option_convert, push_atom!, append_atom!

include("shared.jl")

"""
    AtomicVelocity(atom::Union{AbstractChar,String}, velocity::Vector{Float64})
    AtomicVelocity(x::AtomicSpecies, velocity)
    AtomicVelocity(x::AtomicPosition, velocity)
    AtomicVelocity(atom::Union{AbstractChar,AbstractString})
    AtomicVelocity(x::AtomicSpecies)
    AtomicVelocity(x::AtomicPosition)

Represent each line of the `ATOMIC_VELOCITIES` card in QE's `CP` package.

It is a `mutable struct` and supports _incomplete Initialization_ as in the fourth to
sixth constructors. See the examples below. The `atom` field accepts at most 3 characters
as claimed in QE's documentation.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.CP

julia> AtomicVelocity("H", [0.140374e-04, -0.333683e-04, 0.231834e-04])
AtomicVelocity("H", [1.40374e-5, -3.33683e-5, 2.31834e-5])

julia> x = AtomicVelocity('H');

julia> x.atom
"H"

julia> x.velocity = [0.140374e-04, -0.333683e-04, 0.231834e-04]; x.velocity
3-element Array{Float64,1}:
  1.40374e-5
 -3.33683e-5
  2.31834e-5

julia> AtomicVelocity(AtomicSpecies("H", 1.00794000, "H.pbe-rrkjus_psl.1.0.0.UPF"))
AtomicVelocity("H", #undef)
```
"""
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
AtomicVelocity(atom::AbstractChar, velocity) = AtomicVelocity(string(atom), velocity)
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

"""
    AtomicVelocitiesCard <: Card

Represent the `ATOMIC_VELOCITIES` card in QE's `CP` package. It does not have an "option".
"""
@auto_hash_equals struct AtomicVelocitiesCard <: Card
    data::Vector{AtomicVelocity}
end

"""
    push_atom!(v::AbstractVector{AtomicVelocity}, atoms::AbstractString...)

Push an atom or multiple atoms to a vector of `AtomicVelocity`s.

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

"""
    append_atom!(v::AbstractVector{AtomicVelocity}, atoms::AbstractVector{<:AbstractString})

Append a vector of atoms to a vector of `AtomicVelocity`s.

**Note**: these new `atoms` will result in incomplete `AtomicVelocity`s!

See also: [`append!`](@ref), [`push_atom!`](@ref)
"""
function append_atom!(
    v::AbstractVector{AtomicVelocity},
    atoms::AbstractVector{<:AbstractString},
)
    return append!(v, map(AtomicVelocity, atoms))
end # function append_atom!
"""
    append_atom!(card::AtomicVelocitiesCard, atoms::AbstractVector{<:AbstractString})

Append a vector of atoms to a `AtomicVelocitiesCard`.

**Note**: these new `atoms` will result in incomplete `AtomicVelocity`s!

See also: [`append!`](@ref), [`push_atom!`](@ref)
"""
function append_atom!(card::AtomicVelocitiesCard, atoms::AbstractVector{<:AbstractString})
    append!(card.data, map(AtomicVelocity, atoms))
    return card
end # function append_atom!

@auto_hash_equals struct RefCellParametersCard{T<:Real} <: AbstractCellParametersCard
    option::String
    data::Matrix{T}
    function RefCellParametersCard{T}(option, data) where {T<:Real}
        @assert option ∈ allowed_options(RefCellParametersCard)
        @assert size(data) == (3, 3)
        return new(option, data)
    end
end
RefCellParametersCard(option, data::AbstractMatrix{T}) where {T} =
    RefCellParametersCard{T}(option, data)
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
