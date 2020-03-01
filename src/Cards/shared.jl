using LinearAlgebra: det

using AutoHashEquals: @auto_hash_equals
using Compat: eachrow
using Formatting: sprintf1
using Pseudopotentials: pseudopot_format
using Setfield: get, set, @lens, @set
using StaticArrays: SVector, SMatrix

using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Cards: Card, optionof, allowed_options

import Crystallography
import Pseudopotentials
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Cards

# =============================== AtomicSpecies ============================== #
"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(atom::Union{AbstractChar,AbstractString})
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)
    AtomicSpecies(x::AtomicPosition)

Represent each line of the `ATOMIC_SPECIES` card in QE.

It is a `mutable struct` and supports _incomplete Initialization_ as in the second and
fourth constructors. See the examples below. The `atom` field accepts at most 3 characters
as claimed in QE's documentation.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")

julia> x = AtomicSpecies('S');

julia> x.atom
"S"

julia> x.mass = 32.066; x.mass
32.066

julia> x.pseudopot
ERROR: UndefRefError: access to undefined reference
[...]

julia> x.pseudopot = "S.pz-n-rrkjus_psl.0.1.UPF"; x.pseudopot
"S.pz-n-rrkjus_psl.0.1.UPF"

julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```
"""
@auto_hash_equals mutable struct AtomicSpecies
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
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString})
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        return new(string(atom))
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
# ============================================================================ #

# ============================== AtomicPosition ============================== #
"""
    AtomicPosition(atom::Union{AbstractChar,String}, pos::Vector{Float64}[, if_pos::Vector{Int}])
    AtomicPosition(atom::Union{AbstractChar,AbstractString})
    AtomicPosition(x::AtomicSpecies, pos, if_pos)
    AtomicPosition(x::AtomicSpecies)

Represent each line of the `ATOMIC_POSITIONS` card in QE.

It is a `mutable struct` and supports _incomplete Initialization_ as in the second and
fourth constructors. See the examples below. The `atom` field accepts at most 3 characters
as claimed in QE's documentation.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])

julia> x = AtomicPosition('O');

julia> x.pos
ERROR: UndefRefError: access to undefined reference
[...]

julia> x.pos = [0, 0, 0]
ERROR: TypeError: in setfield!, expected StaticArrays.SArray{Tuple{3},Float64,1,3}, got StaticArrays.SArray{Tuple{3},Int64,1,3}
[...]

julia> x.pos = Float64[0, 0, 0]
3-element Array{Float64,1}:
 0.0
 0.0
 0.0

julia> x.if_pos = [1, 0, 2]
ERROR: TypeError: in setfield!, expected StaticArrays.SArray{Tuple{3},Bool,1,3}, got StaticArrays.SArray{Tuple{3},Int64,1,3}
[...]

julia> x.if_pos = Bool[1, 0, 1]
3-element Array{Bool,1}:
 1
 0
 1

julia> x
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 0, 1])

julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```
"""
@auto_hash_equals mutable struct AtomicPosition
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
    function AtomicPosition(atom::Union{AbstractChar,AbstractString})
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom))
    end
end
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
AtomicPosition(x::AtomicSpecies) = AtomicPosition(x.atom)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicPosition, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)
AtomicSpecies(x::AtomicPosition) = AtomicSpecies(x.atom)

"""
    AtomicPositionsCard <: Card

Represent the `ATOMIC_POSITIONS` card in QE.

# Arguments
- `option::String="alat"`: allowed values are: "alat", "bohr", "angstrom", "crystal", and "crystal_sg".
- `data::AbstractVector{AtomicPosition}`: A vector containing `AtomicPosition`s.
"""
@auto_hash_equals struct AtomicPositionsCard <: Card
    option::String
    data::Vector{AtomicPosition}
    function AtomicPositionsCard(option, data)
        @assert option ∈ allowed_options(AtomicPositionsCard)
        return new(option, data)
    end
end
AtomicPositionsCard(data) = AtomicPositionsCard("alat", data)
function AtomicPositionsCard(option, card::AtomicSpeciesCard)
    return AtomicPositionsCard(option, map(AtomicPosition, card.data))
end # function AtomicPositionsCard
# Introudce mutual constructors since they share the same atoms.
function AtomicSpeciesCard(card::AtomicPositionsCard)
    return AtomicSpeciesCard(map(AtomicSpecies, card.data))
end # function AtomicSpeciesCard

function validate(x::AtomicSpeciesCard, y::AtomicPositionsCard)
    lens = @lens _.data.atom
    @assert(
        isempty(symdiff(map(Base.Fix2(get, lens) ∘ unique, (x, y)))),
        "labels of the atoms are different in `ATOMIC_SPECIES` and `ATOMIC_POSITIONS` card!",
    )
end # function validate
validate(y::AtomicPositionsCard, x::AtomicSpeciesCard) = validate(x, y)

const AtomicSpeciesOrPosition = Union{AtomicSpecies,AtomicPosition}

"""
    push_atom!(v::AbstractVector{AtomicSpecies}, atoms::AbstractString...)
    push_atom!(v::AbstractVector{AtomicPosition}, atoms::AbstractString...)

Push an atom or multiple atoms to a vector of `AtomicSpecies` or `AtomicPosition`s.

**Note**: these new `atoms` will result in incomplete `AtomicSpecies` or `AtomicPosition`s!

See also: [`push!`](@ref), [`append_atom!`](@ref)
"""
function push_atom!(
    v::AbstractVector{T},
    atoms::AbstractString...,
) where {T<:AtomicSpeciesOrPosition}
    return push!(v, map(T, atoms)...)
end # function push_atom!
"""
    push_atom!(card::Union{AtomicSpeciesCard,AtomicPositionsCard}, atoms::AbstractString...)

Push an atom or multiple atoms to a `AtomicSpeciesCard` or `AtomicPositionsCard`.

**Note**: these new `atoms` will result in incomplete `AtomicSpecies` or `AtomicPosition`s!

See also: [`push!`](@ref), [`append_atom!`](@ref)
"""
function push_atom!(
    card::Union{AtomicSpeciesCard,AtomicPositionsCard},
    atoms::AbstractString...,
)
    T = eltype(card.data)
    push!(card.data, map(T, atoms)...)
    return card
end # function push_atom!

"""
    append_atom!(v::AbstractVector{AtomicSpecies}, atoms::AbstractVector{<:AbstractString})
    append_atom!(v::AbstractVector{AtomicPosition}, atoms::AbstractVector{<:AbstractString})

Append a vector of atoms to a vector of `AtomicSpecies` or `AtomicPosition`s.

**Note**: these new `atoms` will result in incomplete `AtomicSpecies` or `AtomicPosition`s!

See also: [`append!`](@ref), [`push_atom!`](@ref)
"""
function append_atom!(
    v::AbstractVector{T},
    atoms::AbstractVector{<:AbstractString},
) where {T<:AtomicSpeciesOrPosition}
    return append!(v, map(T, atoms))
end # function append_atom!
"""
    append_atom!(card::Union{AtomicSpeciesCard,AtomicPositionsCard}, atoms::AbstractVector{<:AbstractString})

Append a vector of atoms to a `AtomicSpeciesCard` or `AtomicPositionsCard`.

**Note**: these new `atoms` will result in incomplete `AtomicSpecies` or `AtomicPosition`s!

See also: [`append!`](@ref), [`push_atom!`](@ref)
"""
function append_atom!(
    card::Union{AtomicSpeciesCard,AtomicPositionsCard},
    atoms::AbstractVector{<:AbstractString},
)
    T = eltype(card.data)
    append!(card.data, map(T, atoms))
    return card
end # function append_atom!
# ============================================================================ #

# ============================== CellParameters ============================== #
"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end

"""
    CellParametersCard{T<:Real} <: AbstractCellParametersCard
    CellParametersCard(option::String, data::Matrix)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
@auto_hash_equals struct CellParametersCard{T<:Real} <: AbstractCellParametersCard
    option::String
    data::SMatrix{3,3,T}
    function CellParametersCard{T}(option, data) where {T<:Real}
        @assert option ∈ allowed_options(CellParametersCard)
        return new(option, data)
    end
end
CellParametersCard(option, data::AbstractMatrix{T}) where {T} =
    CellParametersCard{T}(option, data)
CellParametersCard(data) = CellParametersCard("alat", data)
# ============================================================================ #

# ============================== AtomicForce ============================== #
@auto_hash_equals struct AtomicForce
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
"Represent a general point in the 3D reciprocal space."
abstract type KPoint end

"""
    MonkhorstPackGrid(grid, offsets)

Represent the Monkhorst--Pack grid.

# Arguments
- `grid`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `offsets`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
@auto_hash_equals struct MonkhorstPackGrid
    grid::SVector{3,Int}
    offsets::SVector{3,Bool}
    function MonkhorstPackGrid(grid, offsets)
        # See https://github.com/aiidateam/aiida-quantumespresso/blob/4aef9f9/aiida_quantumespresso/cli/utils/validate.py#L10-L37
        @assert(all(grid .> 0), "`grid` must be positive integers!")
        return new(grid, offsets)
    end
end

"Represent the centre of the Brillouin zone (commonly marked as the Γ point)."
struct GammaPoint <: KPoint end

"""
    SpecialKPoint(coord, weight)

Represent a special point of the 3D Brillouin zone. Each of them has a weight.
"""
@auto_hash_equals struct SpecialKPoint{T<:Real} <: KPoint
    coord::SVector{3,T}
    weight::T
end
SpecialKPoint(coord, weight) = SpecialKPoint{promote_type(eltype(coord), typeof(weight))}(coord, weight)
SpecialKPoint(x, y, z, w) = SpecialKPoint([x, y, z], w)

"""
    struct KPointsCard{<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{<:SpecialKPoint}}} <: Card

Represent the `K_POINTS` card in QE.

# Arguments
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
- `data::Union{MonkhorstPackGrid,GammaPoint,AbstractVector{<:SpecialKPoint}}`: A Γ point, a Monkhorst--Pack grid or a vector containing `SpecialKPoint`s.
"""
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

# function meshgrid(reciprocal::AbstractMatrix, mp::MonkhorstPackGrid, crystal::Bool = true)
#     nk, sk = mp.grid, mp.offsets
#     sk = map(x -> isone(x) ? 1 // 2 : 0, sk)
#     mesh = Iterators.product(
#         (0:nk[1]-1) // nk[1] .+ sk[1],
#         (0:nk[2]-1) // nk[2] .+ sk[2],
#         (0:nk[3]-1) // nk[3] .+ sk[3],
#     )
#     if crystal
#         return collect(mesh)
#     else
#         a1, a2, a3 = reciprocal[1, :], reciprocal[2, :], reciprocal[3, :]
#         return collect(x[1] .* a1 + x[2] .* a2 + x[3] .* a3 for x in mesh)
#     end
# end # function meshgrid
# function meshgrid(card::CellParametersCard, mp::MonkhorstPackGrid)
#     option = optionof(card)
#     if option == "alat"
#         return meshgrid(reciprocalof(card.data), mp, true)
#     else  # option ∈ ("bohr", "angstrom")
#         return meshgrid(reciprocalof(card.data), mp, false)
#     end
# end # function meshgrid
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
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

It will throw an error if the information is not enough to calculate the volume.
"""
function Crystallography.cellvolume(card::AbstractCellParametersCard)
    BOHR_TO_ANGSTROM = 0.529177210903
    option = optionof(card)
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
    option_convert(new_option::AbstractString, card::AbstractCellParametersCard)

Convert the option of an `AbstractCellParametersCard` from "bohr" to "angstrom", or its reverse.

It does not support conversion between "alat" and the rests.
"""
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
    newline = '\n',
)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES" *
    newline *
    join((indent * to_qe(x; delim = delim, numfmt = numfmt) for x in card.data), newline)
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
    newline = '\n',
    verbose::Bool = false,
)
    return "ATOMIC_POSITIONS { $(optionof(card)) }" *
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
    return "CELL_PARAMETERS { $(optionof(card)) }" *
    newline *
    join(
        (
            indent * join(map(x -> sprintf1(numfmt, x), row), delim)
            for row in eachrow(card.data)
        ),
        newline,
    )
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
    if optionof(card) in ("gamma", "automatic")
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
