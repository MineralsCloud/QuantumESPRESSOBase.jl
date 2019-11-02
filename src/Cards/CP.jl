module CP

using parameters: @with_kw
using QuantumESPRESSOBase
using ..Cards

export AtomicVelocity,
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
@with_kw struct RefCellParametersCard{A<:AbstractVector,B<:AbstractMatrix} <: Card
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

"""
    option(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `option(x)` is suggested.
"""
Cards.option(::AtomicVelocitiesCard) = "a.u"
Cards.option(::AtomicForcesCard) = nothing

"""
    allowed_options(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards, QuantumESPRESSOBase.Cards.PWscf

julia> allowed_options(AtomicVelocitiesCard)
("a.u",)

julia> allowed_options(RefCellParametersCard)
("bohr", "angstrom")
```
"""
Cards.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Cards.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

const ANGSTROM_TO_BOHR = 1 / 0.529177210903

"""
    ref_cell_volume(card::RefCellParametersCard)

Return the referece cell volume according to the `RefCellParametersCard`'s parameters, in atomic unit.
"""
function ref_cell_volume(card::RefCellParametersCard)
    @match option(card) begin
        "bohr" => det(card.data)
        "angstrom" => det(card.data) * ANGSTROM_TO_BOHR^3
    end
end # function ref_cell_volume

end
