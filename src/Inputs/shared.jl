using AutoHashEquals: @auto_hash_equals
using Compat: eachrow
using Crystallography: Lattice, Cell
using Formatting: sprintf1
using Pseudopotentials: pseudopot_format
using Setfield: get, set, @lens, @set
using StaticArrays: SVector, SMatrix, FieldVector

using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Inputs: Card, getoption, allowed_options

import Crystallography
import Pseudopotentials
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Inputs

# =============================== AtomicSpecies ============================== #
"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)

Represent each line of the `ATOMIC_SPECIES` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")

julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```
"""
struct AtomicSpecies
    "Label of the atom. Max total length cannot exceed 3 characters."
    atom::String
    """
    Mass of the atomic species in atomic unit.

    Used only when performing molecular dynamics (MD) run
    or structural optimization runs using damped MD.
    Not actually used in all other cases (but stored
    in data files, so phonon calculations will use
    these values unless other values are provided).
    """
    mass::Float64
    """
    File containing pseudopotential for this species.

    See also: [`pseudopot_format`](@ref)
    """
    pseudopot::String
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, pseudopot)
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        return new(string(atom), mass, pseudopot)
    end
end

"""
    pseudopot_format(data::AtomicSpecies)

Return the pseudopotential format of the `AtomicSpecies`.
"""
Pseudopotentials.pseudopot_format(data::AtomicSpecies) = pseudopot_format(data.pseudopot)

"""
    AtomicSpeciesCard <: Card

Represent the `ATOMIC_SPECIES` card in QE. It does not have an "option".
"""
struct AtomicSpeciesCard <: Card
    data::Vector{AtomicSpecies}
end
AtomicSpeciesCard(cell::Cell) = AtomicSpeciesCard(map(AtomicSpecies ∘ string, cell.numbers))
# ============================================================================ #

# ============================== AtomicPosition ============================== #
"""
    AtomicPosition(atom::Union{AbstractChar,String}, pos::Vector{Float64}[, if_pos::Vector{Int}])
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Represent each line of the `ATOMIC_POSITIONS` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])

julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```
"""
struct AtomicPosition
    "Label of the atom as specified in `AtomicSpecies`."
    atom::String
    "Atomic positions. A three-element vector of floats."
    pos::SVector{3,Float64}
    """
    Component `i` of the force for this atom is multiplied by `if_pos(i)`,
    which must be either `0` or `1`.  Used to keep selected atoms and/or
    selected components fixed in MD dynamics or structural optimization run.

    With `crystal_sg` atomic coordinates the constraints are copied in all equivalent
    atoms.
    """
    if_pos::SVector{3,Bool}
    function AtomicPosition(atom::Union{AbstractChar,AbstractString}, pos, if_pos)
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom), pos, if_pos)
    end
end
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicPosition, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)

"""
    AtomicPositionsCard <: Card

Represent the `ATOMIC_POSITIONS` card in QE.

# Arguments
- `data::AbstractVector{AtomicPosition}`: A vector containing `AtomicPosition`s.
- `option::String="alat"`: allowed values are: "alat", "bohr", "angstrom", "crystal", and "crystal_sg".
"""
struct AtomicPositionsCard <: Card
    data::Vector{AtomicPosition}
    option::String
    function AtomicPositionsCard(data, option = "alat")
        @assert option ∈ allowed_options(AtomicPositionsCard)
        return new(data, option)
    end
end
AtomicPositionsCard(cell::Cell, option) = AtomicPositionsCard([AtomicPosition(string(atom), pos) for (atom, pos) in zip(cell.numbers, cell.positions)], option)
# Introudce mutual constructors since they share the same atoms.

function validate(x::AtomicSpeciesCard, y::AtomicPositionsCard)
    lens = @lens _.data.atom
    @assert(
        isempty(symdiff(map(Base.Fix2(get, lens) ∘ unique, (x, y)))),
        "labels of the atoms are different in `ATOMIC_SPECIES` and `ATOMIC_POSITIONS` card!",
    )
end # function validate
validate(y::AtomicPositionsCard, x::AtomicSpeciesCard) = validate(x, y)

# ============================== CellParameters ============================== #
"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end

"""
    CellParametersCard{T<:Real} <: AbstractCellParametersCard
    CellParametersCard(data::AbstractMatrix, option::String)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
struct CellParametersCard{T<:Real} <: AbstractCellParametersCard
    data::SMatrix{3,3,T}
    option::String
    function CellParametersCard{T}(data, option = "alat") where {T<:Real}
        @assert option ∈ allowed_options(CellParametersCard)
        return new(data, option)
    end
end
CellParametersCard(data::AbstractMatrix{T}, option = "alat") where {T} =
    CellParametersCard{T}(data, option)
CellParametersCard(lattice::Lattice{T}, option) where {T} = CellParametersCard(convert(Matrix{T}, lattice), option)
CellParametersCard(cell::Cell, option) = CellParametersCard(cell.lattice, option)
# ============================================================================ #

# ============================== AtomicForce ============================== #
struct AtomicForce
    atom::String
    force::SVector{3,Float64}
    function AtomicForce(atom::Union{AbstractChar,AbstractString}, force)
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom), force)
    end
end

struct AtomicForcesCard <: Card
    data::Vector{AtomicForce}
end
# ============================================================================ #

# ============================== KPointsCard ============================== #
"""
    MonkhorstPackGrid(grid, offsets)

Represent the Monkhorst--Pack grid.

# Arguments
- `grid`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `offsets`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    grid::SVector{3,UInt}
    offsets::SVector{3,Bool}
end

"Represent the centre of the Brillouin zone (commonly marked as the Γ point)."
struct GammaPoint end

"""
    SpecialKPoint(coord, weight)

Represent a special point of the 3D Brillouin zone. Each of them has a weight.
"""
struct SpecialKPoint <: FieldVector{4,Float64}
    x::Float64
    y::Float64
    z::Float64
    w::Float64
end

"""
    struct KPointsCard{<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}} <: Card

Represent the `K_POINTS` card in QE.

# Arguments
- `data::Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}`: A Γ point, a Monkhorst--Pack grid or a vector containing `SpecialKPoint`s.
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
"""
struct KPointsCard{
    A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}},
} <: Card
    data::A
    option::String
    function KPointsCard{A}(
        data,
        option,
    ) where {A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}}
        @assert option ∈ allowed_options(KPointsCard)
        @assert if option == "automatic"
            typeof(data) <: MonkhorstPackGrid
        elseif option == "gamma"
            typeof(data) <: GammaPoint
        else  # option ∈ ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            eltype(data) <: SpecialKPoint
        end
        return new(data, option)
    end
end
KPointsCard(data::A, option) where {A} = KPointsCard{A}(data, option)
KPointsCard(data::AbstractVector{SpecialKPoint}) = KPointsCard(data, "tpiba")
KPointsCard(data::GammaPoint) = KPointsCard(data, "gamma")
KPointsCard(data::MonkhorstPackGrid) = KPointsCard(data, "automatic")

Inputs.getoption(::AtomicSpeciesCard) = nothing

Inputs.allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
Inputs.allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
Inputs.allowed_options(::Type{<:KPointsCard}) = (
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
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

It will throw an error if the information is not enough to calculate the volume.
"""
function Crystallography.cellvolume(card::AbstractCellParametersCard)
    BOHR_TO_ANGSTROM = 0.529177210903
    option = getoption(card)
    if option == "bohr"
        abs(det(card.data))
    elseif option == "angstrom"
        abs(det(card.data)) / BOHR_TO_ANGSTROM^3
    elseif option == "alat"
        error("Information not enough! The `celldm[1]` parameter is unknown!")
    else
        error("Option $option is unknown!")
    end
end # function Crystallography.cellvolume

"""
    optconvert(new_option::AbstractString, card::AbstractCellParametersCard)

Convert the option of an `AbstractCellParametersCard` from "bohr" to "angstrom", or its reverse.

It does not support conversion between "alat" and the rests.
"""
function optconvert(new_option::AbstractString, card::AbstractCellParametersCard)
    BOHR_TO_ANGSTROM = 0.529177210903
    old_option = getoption(card)
    new_option == old_option && return card  # No conversion is needed
    pair = old_option => new_option
    factor = if pair == ("bohr" => "angstrom")
        BOHR_TO_ANGSTROM
    elseif pair == ("angstrom" => "bohr")
        1 / BOHR_TO_ANGSTROM
    else
        error("Unknown option pair ($pair) given!")
    end
    return typeof(card)(card.data .* factor, new_option)
end # function optconvert

Inputs.entryname(::Type{ControlNamelist}) = :control
Inputs.entryname(::Type{SystemNamelist}) = :system
Inputs.entryname(::Type{ElectronsNamelist}) = :electrons
Inputs.entryname(::Type{IonsNamelist}) = :ions
Inputs.entryname(::Type{CellNamelist}) = :cell
Inputs.entryname(::Type{AtomicSpeciesCard}) = :atomic_species
Inputs.entryname(::Type{AtomicPositionsCard}) = :atomic_positions
Inputs.entryname(::Type{<:CellParametersCard}) = :cell_parameters
Inputs.entryname(::Type{<:KPointsCard}) = :k_points

Inputs.titleof(::Type{ControlNamelist}) = "CONTROL"
Inputs.titleof(::Type{SystemNamelist}) = "SYSTEM"
Inputs.titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
Inputs.titleof(::Type{IonsNamelist}) = "IONS"
Inputs.titleof(::Type{CellNamelist}) = "CELL"
Inputs.titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
Inputs.titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
Inputs.titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"
Inputs.titleof(::Type{<:KPointsCard}) = "K_POINTS"

"""
    Crystallography.BravaisLattice(nml::SystemNamelist)

Return a `BravaisLattice` from a `SystemNamelist`.
"""
Crystallography.BravaisLattice(nml::SystemNamelist) = BravaisLattice{nml.ibrav}()

"""
    Crystallography.Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
function Crystallography.Lattice(nml::SystemNamelist)
    b = BravaisLattice(nml)
    return Lattice(b, CellParameters(nml.celldm...))
end # function Crystallography.Lattice

"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
Crystallography.cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))

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
    newline = '\n',
)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES" *
    newline *
    join((indent * to_qe(x; delim = delim, numfmt = numfmt) for x in unique(card.data)), newline)
end
function QuantumESPRESSOBase.to_qe(
    data::AtomicPosition;
    delim = ' ',
    numfmt = "%14.9f",
    verbose::Bool = false,
)
    v = verbose ? [data.pos; data.if_pos] : data.pos
    return join([sprintf1("%3s", data.atom); map(x -> sprintf1(numfmt, x), v)], delim)  # FIXME: `numfmt` on `Bool` will give wrong result
end
function QuantumESPRESSOBase.to_qe(
    card::AtomicPositionsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
    verbose::Bool = false,
)
    return "ATOMIC_POSITIONS { $(getoption(card)) }" *
    newline *
    join(
        (
            indent * to_qe(x; delim = delim, numfmt = numfmt, verbose = verbose)
            for x in card.data
        ),
        newline,
    )
end
function QuantumESPRESSOBase.to_qe(
    card::CellParametersCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    it = (indent * join((sprintf1(numfmt, x) for x in row), delim) for row in eachrow(card.data))
    return "CELL_PARAMETERS { $(getoption(card)) }" * newline * join(it, newline)
end
QuantumESPRESSOBase.to_qe(data::GammaPoint) = ""
function QuantumESPRESSOBase.to_qe(data::MonkhorstPackGrid; delim = ' ', numfmt = "%5d")
    return join(
        map(x -> sprintf1(numfmt, x), [getfield(data, 1); getfield(data, 2)]),
        delim,
    )
end
function QuantumESPRESSOBase.to_qe(data::SpecialKPoint; delim = ' ', numfmt = "%14.9f")
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
    newline = '\n',
)
    content = "K_POINTS { $(card.option) }" * newline
    if getoption(card) in ("gamma", "automatic")
        content *= indent * to_qe(card.data)
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= string(length(card.data), newline)
        content *= join(
            (indent * to_qe(x; delim = delim, numfmt = numfmt) for x in card.data),
            newline,
        )
    end
    return content
end

function Base.setproperty!(value::AtomicSpecies, name::Symbol, x)
    if name == :atom
        @assert(length(x) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        x = string(x)  # An `if` statement is more expensive than directly setting a string
    end
    setfield!(value, name, x)
end # function Base.setproperty!
function Base.setproperty!(value::AtomicPosition, name::Symbol, x)
    x = if name == :atom
        @assert(length(x) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        x = string(x)  # Make sure it is a `String`
    elseif name ∈ (:pos, :if_pos) && x isa AbstractVector
        SVector{3}(x)
    end
    setfield!(value, name, x)
end # function Base.setproperty!
