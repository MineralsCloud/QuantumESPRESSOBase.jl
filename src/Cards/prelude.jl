abstract type Card <: InputEntry end
abstract type AbstractCellParametersCard <: Card end

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies
    atom::String
    mass::Float64
    potential::String
end

abstract type PseudopotentialFormat end
struct VanderbiltUltraSoft <: PseudopotentialFormat end
struct AndreaDalCorso <: PseudopotentialFormat end
struct OldNormConserving <: PseudopotentialFormat end

"""
    potential_format(data::AtomicSpecies)::String

Return the pseudopotential format.

The pseudopotential file is assumed to be in the new UPF format.
If it doesn't work, the pseudopotential format is determined by
the file name:
- "*.vdb or *.van": Vanderbilt US pseudopotential code
- "*.RRKJ3": Andrea Dal Corso's code (old format)
- none of the above: old PWscf norm-conserving format
"""
function potential_format(data::AtomicSpecies)::PseudopotentialFormat
    @match lowercase(splitext(data.pseudo)[2]) begin
        ".vdb" || ".van" => VanderbiltUltraSoft()
        ".rrkj3" => AndreaDalCorso()
        _ => OldNormConserving()
    end
end

struct AtomicSpeciesCard{T<:AbstractVector{<:AtomicSpecies}} <: Card
    data::T
end
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{A<:AbstractVector{<:Real},B<:AbstractVector{<:Integer}}
    atom::String
    pos::A
    if_pos::B = [1, 1, 1]
    @assert(length(pos) == 3, "`pos` is not of length 3, but $(length(pos))!",)
    @assert(length(if_pos) == 3, "`if_pos` is not of length 3, but $(length(if_pos))!",)
    @assert(all(x ∈ (0, 1) for x in if_pos), "`if_pos` must be either 0 or 1!")
end
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, [1, 1, 1])

@with_kw struct AtomicPositionsCard{A<:AbstractVector{<:AtomicPosition}} <: Card
    option::String = "alat"
    data::A
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
@with_kw struct CellParametersCard{A<:AbstractMatrix{<:Real}} <: AbstractCellParametersCard
    option::String = "alat"
    data::A
    @assert(option ∈ allowed_options(CellParametersCard))
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

"""
    optionof(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `optionof(x)` is suggested.
"""
optionof(card::Card) = getfield(card, :option)
optionof(::AtomicSpeciesCard) = nothing

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

const BOHR_TO_ANGSTROM = 0.529177210903

"""
    cell_volume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.
"""
function cell_volume(card::AbstractCellParametersCard)
    @match optionof(card) begin
        "bohr" => det(card.data)
        "angstrom" => det(card.data) / BOHR_TO_ANGSTROM^3
        "alat" => error("Information not enough! The `celldm[1]` parameter is unknown!")
    end
end # function cell_volume

function option_convert(new_option::AbstractString, card::AbstractCellParametersCard)
    old_option = optionof(card)
    factor = @match (old_option => new_option) begin
        ("bohr" => "angstrom") => BOHR_TO_ANGSTROM
        ("angstrom" => "bohr") => 1 / BOHR_TO_ANGSTROM
        _ => error("Unknown option pair ($old_option => $new_option) given!")
    end
    return typeof(card)(new_option, card.data .* factor)
end # function option_convert
