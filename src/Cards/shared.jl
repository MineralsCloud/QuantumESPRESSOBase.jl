using LinearAlgebra: det

using AutoHashEquals: @auto_hash_equals
using Compat: eachrow
using Formatting: sprintf1
using Setfield: @lens, get

using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Cards: Card, optionof, allowed_options

import QuantumESPRESSOBase
import QuantumESPRESSOBase.Cards

# =============================== AtomicSpecies ============================== #
@auto_hash_equals struct AtomicSpecies
    atom::String
    mass::Float64
    pseudopot::String
    function AtomicSpecies(atom, mass, pseudopot)
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        return new(atom, mass, pseudopot)
    end
end
AtomicSpecies(atom::AbstractChar, mass, pseudopot) =
    AtomicSpecies(string(atom), mass, pseudopot)

abstract type PseudopotentialFormat end
"""
    UnifiedPseudopotentialFormat <: PseudopotentialFormat

A singleton representing the new UPF format.

If it doesn't work, the pseudopotential format is determined by
the file name.
"""
struct UnifiedPseudopotentialFormat <: PseudopotentialFormat end
"""
    VanderbiltUltraSoft <: PseudopotentialFormat

A singleton representing the Vanderbilt US pseudopotential code.
"""
struct VanderbiltUltraSoft <: PseudopotentialFormat end
"""
    AndreaDalCorso <: PseudopotentialFormat

A singleton representing the Andrea Dal Corso's code (old format).
"""
struct AndreaDalCorso <: PseudopotentialFormat end
"""
    OldNormConserving <: PseudopotentialFormat

A singleton representing the old PWscf norm-conserving format.
"""
struct OldNormConserving <: PseudopotentialFormat end

"""
    pseudopot_format(data::AtomicSpecies)::String

Return the pseudopotential format.

The pseudopotential file is assumed to be in the new UPF format.
If it doesn't work, the pseudopotential format is determined by
the file name:
- "*.vdb or *.van": Vanderbilt US pseudopotential code
- "*.RRKJ3": Andrea Dal Corso's code (old format)
- none of the above: old PWscf norm-conserving format
"""
function pseudopot_format(data::AtomicSpecies)::PseudopotentialFormat
    ext = uppercase(splitext(data.pseudopot)[2])
    return if ext == ".UPF"
        UnifiedPseudopotentialFormat()
    elseif ext ∈ (".VDB", ".VAN")
        VanderbiltUltraSoft()
    elseif ext == ".RRKJ3"
        AndreaDalCorso()
    else
        OldNormConserving()
    end
end

struct AtomicSpeciesCard{T<:AbstractVector{<:AtomicSpecies}} <: Card
    data::T
end
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@auto_hash_equals struct AtomicPosition{
    A<:AbstractVector{<:Real},
    B<:AbstractVector{<:Integer},
}
    atom::String
    pos::A
    if_pos::B
    function AtomicPosition{A,B}(
        atom,
        pos,
        if_pos,
    ) where {A<:AbstractVector{<:Real},B<:AbstractVector{<:Integer}}
        @assert length(pos) == length(if_pos) == 3
        @assert(
            all(iszero(x) || isone(x) for x in if_pos),
            "`if_pos` elements must be 0 or 1!"
        )
        return new(atom, pos, if_pos)
    end
end
AtomicPosition(atom, pos::A, if_pos::B) where {A,B} = AtomicPosition{A,B}(atom, pos, if_pos)
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, ones(Int, 3))
AtomicPosition(x::AbstractChar, pos, if_pos) = AtomicPosition(string(x), pos, if_pos)
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)

@auto_hash_equals struct AtomicPositionsCard{A<:AbstractVector{<:AtomicPosition}} <: Card
    option::String
    data::A
    function AtomicPositionsCard{A}(
        option,
        data,
    ) where {A<:AbstractVector{<:AtomicPosition}}
        @assert option ∈ allowed_options(AtomicPositionsCard)
        return new(option, data)
    end
end
AtomicPositionsCard(option, data::A) where {A} = AtomicPositionsCard{A}(option, data)
AtomicPositionsCard(data) = AtomicPositionsCard("alat", data)

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
abstract type AbstractCellParametersCard <: Card end
@auto_hash_equals struct CellParametersCard{A<:AbstractMatrix{<:Real}} <:
                         AbstractCellParametersCard
    option::String
    data::A
    function CellParametersCard{A}(option, data) where {A<:AbstractMatrix{<:Real}}
        @assert option ∈ allowed_options(CellParametersCard)
        @assert size(data) == (3, 3)
        return new(option, data)
    end
end
CellParametersCard(option, data::A) where {A} = CellParametersCard{A}(option, data)
CellParametersCard(data) = CellParametersCard("alat", data)
# ============================================================================ #

# ============================== AtomicForce ============================== #
@auto_hash_equals struct AtomicForce{A<:AbstractVector{<:Real}}
    atom::String
    force::A
    function AtomicForce{A}(atom, force) where {A<:AbstractVector{<:Real}}
        @assert length(force) == 3
        return new(atom, force)
    end
end
AtomicForce(atom, force::A) where {A} = AtomicForce{A}(atom, force)

struct AtomicForcesCard{T<:AbstractVector{<:AtomicForce}} <: Card
    data::T
end
# ============================================================================ #

# ============================== KPointsCard ============================== #
abstract type KPoint end

@auto_hash_equals struct MonkhorstPackGrid{
    A<:AbstractVector{<:Integer},
    B<:AbstractVector{<:Integer},
}
    grid::A
    offsets::B
    function MonkhorstPackGrid{A,B}(grid, offsets) where {A,B}
        @assert length(grid) == length(offsets) == 3
        # See https://github.com/aiidateam/aiida-quantumespresso/blob/4aef9f9/aiida_quantumespresso/cli/utils/validate.py#L10-L37
        @assert(all(grid .> 0), "`grid` must be positive integers!")
        @assert(
            all(iszero(x) || isone(x) for x in offsets),
            "`offsets` elements must be 0 or 1!"
        )
        return new(grid, offsets)
    end
end
MonkhorstPackGrid(grid::A, offsets::B) where {A,B} = MonkhorstPackGrid{A,B}(grid, offsets)

struct GammaPoint <: KPoint end

@auto_hash_equals struct SpecialKPoint{A<:AbstractVector{<:Real},B<:Real} <: KPoint
    coord::A
    weight::B
    function SpecialKPoint{A,B}(coord, weight) where {A<:AbstractVector{<:Real},B<:Real}
        @assert length(coord) == 3
        return new(coord, weight)
    end
end
SpecialKPoint(coord::A, weight::B) where {A,B} = SpecialKPoint{A,B}(coord, weight)
SpecialKPoint(x, y, z, w) = SpecialKPoint([x, y, z], w)

@auto_hash_equals struct KPointsCard{
    A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{<:SpecialKPoint}},
} <: Card
    option::String
    data::A
    function KPointsCard{A}(
        option,
        data,
    ) where {A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{<:SpecialKPoint}}}
        @assert option ∈ allowed_options(KPointsCard)
        @assert if option == "automatic"
            typeof(data) <: MonkhorstPackGrid
        elseif option == "gamma"
            typeof(data) <: GammaPoint
        else  # option ∈ ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            eltype(data) <: SpecialKPoint
        end
        return new(option, data)
    end
end
KPointsCard(option, data::A) where {A} = KPointsCard{A}(option, data)
KPointsCard(data) = KPointsCard("tpiba", data)
function KPointsCard(option::AbstractString, data::AbstractMatrix{<:Real})
    @assert(size(data, 2) == 4, "The size of `data` is not `(N, 4)`, but $(size(data))!")
    return KPointsCard(option, [SpecialKPoint(x...) for x in eachrow(data)])
end
# ============================================================================ #

Cards.optionof(::AtomicSpeciesCard) = nothing

Cards.allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
Cards.allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
Cards.allowed_options(::Type{<:KPointsCard}) = (
    "tpiba",
    "automatic",
    "crystal",
    "gamma",
    "tpiba_b",
    "crystal_b",
    "tpiba_c",
    "crystal_c",
)

"""
    cell_volume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.
"""
function QuantumESPRESSOBase.cell_volume(card::AbstractCellParametersCard)
    BOHR_TO_ANGSTROM = 0.529177210903
    option = optionof(card)
    if option == "bohr"
        det(card.data)
    elseif option == "angstrom"
        det(card.data) / BOHR_TO_ANGSTROM^3
    elseif option == "alat"
        error("Information not enough! The `celldm[1]` parameter is unknown!")
    else
        error("Option $option is unknown!")
    end
end # function QuantumESPRESSOBase.cell_volume

function option_convert(new_option::AbstractString, card::AbstractCellParametersCard)
    BOHR_TO_ANGSTROM = 0.529177210903
    old_option = optionof(card)
    new_option == old_option && return card  # No conversion is needed
    pair = old_option => new_option
    factor = if pair == ("bohr" => "angstrom")
        BOHR_TO_ANGSTROM
    elseif pair == ("angstrom" => "bohr")
        1 / BOHR_TO_ANGSTROM
    else
        error("Unknown option pair ($pair) given!")
    end
    return typeof(card)(new_option, card.data .* factor)
end # function option_convert

QuantumESPRESSOBase.asfieldname(::Type{<:AtomicSpeciesCard}) = :atomic_species
QuantumESPRESSOBase.asfieldname(::Type{<:AtomicPositionsCard}) = :atomic_positions
QuantumESPRESSOBase.asfieldname(::Type{<:CellParametersCard}) = :cell_parameters
QuantumESPRESSOBase.asfieldname(::Type{<:KPointsCard}) = :k_points

QuantumESPRESSOBase.titleof(::Type{<:AtomicSpeciesCard}) = "ATOMIC_SPECIES"
QuantumESPRESSOBase.titleof(::Type{<:AtomicPositionsCard}) = "ATOMIC_POSITIONS"
QuantumESPRESSOBase.titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"
QuantumESPRESSOBase.titleof(::Type{<:KPointsCard}) = "K_POINTS"

function QuantumESPRESSOBase.to_qe(data::AtomicSpecies; delim = ' ', numfmt = "%14.9f")
    return join(
        (sprintf1("%3s", data.atom), sprintf1(numfmt, data.mass), data.pseudopot),
        delim,
    )
end
function QuantumESPRESSOBase.to_qe(
    card::AtomicSpeciesCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%20.10f",
)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES\n" * join(
        (indent * to_qe(x; delim = delim, numfmt = numfmt) for x in card.data),
        "\n",
    )
end
function QuantumESPRESSOBase.to_qe(
    data::AtomicPosition;
    delim = ' ',
    numfmt = "%14.9f",
    verbose::Bool = false,
)
    v = verbose ? [data.pos; data.if_pos] : data.pos
    return join([sprintf1("%3s", data.atom); map(x -> sprintf1(numfmt, x), v)], delim)
end
function QuantumESPRESSOBase.to_qe(
    card::AtomicPositionsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    verbose::Bool = false,
)
    return "ATOMIC_POSITIONS { $(optionof(card)) }\n" * join(
        (
            indent * to_qe(x; delim = delim, numfmt = numfmt, verbose = verbose)
            for x in card.data
        ),
        "\n",
    )
end
function QuantumESPRESSOBase.to_qe(
    card::CellParametersCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
)
    return "CELL_PARAMETERS { $(optionof(card)) }\n" * join(
        (
            indent * join(map(x -> sprintf1(numfmt, x), row), delim)
            for row in eachrow(card.data)
        ),
        "\n",
    )
end
QuantumESPRESSOBase.to_qe(data::GammaPoint) = ""
function QuantumESPRESSOBase.to_qe(
    data::Union{MonkhorstPackGrid,SpecialKPoint};
    delim = ' ',
    numfmt = "%14.9f",
)
    return join(
        map(x -> sprintf1(numfmt, x), [getfield(data, 1); getfield(data, 2)]),
        delim,
    )
end
function QuantumESPRESSOBase.to_qe(
    card::KPointsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
)
    content = "K_POINTS { $(card.option) }\n"
    if optionof(card) in ("gamma", "automatic")
        content *= indent * to_qe(card.data) * "\n"
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= "$(length(card.data))\n"
        for x in card.data
            content *= indent * to_qe(x; delim = delim, numfmt = numfmt) * "\n"
        end
    end
    return content
end