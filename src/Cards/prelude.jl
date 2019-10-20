abstract type Card <: InputEntry end

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies{A<:AbstractString,B<:Real,C<:AbstractString}
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
        ".vdb" || ".van" => "Vanderbilt US pseudopotential code"
        ".rrkj3" => "Andrea Dal Corso's code (old format)"
        _ => "old PWscf norm-conserving format"
    end
end

struct AtomicSpeciesCard{T<:AbstractVector{<:AtomicSpecies}} <: Card
    data::T
end
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{
    A<:AbstractString,
    B<:AbstractVector{<:Real},
    C<:AbstractVector{Int},
}
    atom::A
    pos::B
    if_pos::C = [1, 1, 1]
    @assert(
        length(pos) == 3,
        "`pos` must be a three-element-vector! However it is of length $(length(pos))!",
    )
    @assert(
        length(if_pos) == 3,
        "`if_pos` must be a three-element-vector! However it is of length $(length(if_pos))!",
    )
    @assert(all(x ∈ (0, 1) for x in if_pos), "`if_pos` must be either 0 or 1!")
end
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, [1, 1, 1])

@with_kw struct AtomicPositionsCard{
    A<:AbstractString,
    B<:AbstractVector{<:AtomicPosition},
} <: Card
    option::A = "alat"
    data::B
    @assert(option ∈ allowed_options(AtomicPositionsCard))
end

function validate(x::AtomicSpeciesCard, y::AtomicPositionsCard)
    lens = @lens _.data.atom
    @assert(
        isempty(symdiff(map(Base.Fix2(get, lens) ∘ unique, (x, y)))),
        "labels of the atoms are different in `ATOMIC_SPECIES` and `ATOMIC_POSITIONS` card!",
    )
end # function validate
validate(y::AtomicPositionsCard, x::AtomicSpeciesCard) = validate(x, y)
# ============================================================================ #

# ============================== CellParameters ============================== #
@with_kw struct CellParametersCard{A<:AbstractString,B<:AbstractMatrix} <: Card
    option::A = "alat"
    data::B
    @assert(option ∈ allowed_options(CellParametersCard))
    @assert(size(data) == (3, 3))
end
# ============================================================================ #

"""
    option(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `option(x)` is suggested.
"""
option(card::Card) = getfield(card, :option)
option(card::AtomicSpeciesCard) = nothing

"""
    allowed_options(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards, QuantumESPRESSOBase.Cards.PWscf

julia> allowed_options(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> allowed_options(CellParametersCard)
("alat", "bohr", "angstrom")

julia> allowed_options(KPointsCard)
("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
allowed_options(::Type{<:Card}) = nothing
allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
