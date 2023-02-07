using Crystallography: Cell, ReciprocalPoint, MonkhorstPackGrid
using Functors: fmap
using StaticArrays: SVector, SMatrix

import ..Inputs: optionpool

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    KPointsCard,
    KMeshCard,
    GammaPointCard,
    SpecialPointsCard
export optconvert, optionpool, eachatom, getpotentials

"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)

Represent each line of the `ATOMIC_SPECIES` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```julia
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

    See also: [`pseudoformat`](@ref)
    """
    pseudopot::String
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, pseudopot)
        @assert length(atom) <= 3 "`atom` can have at most 3 characters!"
        return new(string(atom), mass, pseudopot)
    end
end

"""
    pseudoformat(data::AtomicSpecies)

Return the pseudopotential format of the `AtomicSpecies`.
"""
# pseudoformat(data::AtomicSpecies) = pseudoformat(data.pseudopot)

"""
    AtomicSpeciesCard <: Card

Represent the `ATOMIC_SPECIES` card in QE. It does not have an "option".
"""
@struct_hash_equal struct AtomicSpeciesCard <: Card
    data::Vector{AtomicSpecies}
end

"""
    getpotentials(card::AtomicSpeciesCard)

Get the pseudopotential names from an `AtomicSpeciesCard`.
"""
function getpotentials(card::AtomicSpeciesCard)
    return map(card.data) do atomic_species
        atomic_species.pseudopot
    end
end

"""
    AtomicPosition(atom::Union{AbstractChar,String}, pos::Vector{Float64}[, if_pos::Vector{Int}])
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Represent each line of the `ATOMIC_POSITIONS` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```julia
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
        @assert length(atom) <= 3 "`atom` can have at most 3 characters!"
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
@struct_hash_equal struct AtomicPositionsCard <: Card
    data::Vector{AtomicPosition}
    option::String
    function AtomicPositionsCard(data, option="alat")
        @assert option in optionpool(AtomicPositionsCard)
        return new(data, option)
    end
end
AtomicPositionsCard(cell::Cell, option) = AtomicPositionsCard(
    map(cell.types, eachcol(cell.positions)) do atom, position
        AtomicPosition(string(atom), position)
    end,
    option,
)

"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end

"""
    CellParametersCard{T<:Real} <: AbstractCellParametersCard
    CellParametersCard(data::AbstractMatrix, option::String)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
struct CellParametersCard <: AbstractCellParametersCard
    data::SMatrix{3,3,Float64}
    option::String
    function CellParametersCard(data, option="alat")
        @assert option in optionpool(CellParametersCard)
        return new(data, option)
    end
end
CellParametersCard(lattice::Lattice{<:Real}, option) =
    CellParametersCard(transpose(lattice.data), option)
CellParametersCard(lattice::Lattice{<:Length}) =
    CellParametersCard(fmap(x -> ustrip(u"bohr", x), lattice), "bohr")
CellParametersCard(cell::Cell, option) = CellParametersCard(transpose(cell.lattice), option)

struct AtomicForce
    atom::String
    force::SVector{3,Float64}
    function AtomicForce(atom::Union{AbstractChar,AbstractString}, force)
        @assert length(atom) <= 3 "`atom` can have at most 3 characters!"
        return new(string(atom), force)
    end
end

@struct_hash_equal struct AtomicForcesCard <: Card
    data::Vector{AtomicForce}
end

# See https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1008-L1032 & https://github.com/JuliaLang/julia/blob/de3a70a/base/io.jl#L971-L1054
struct EachAtom{T}
    card::T
end

eachatom(card::Union{AtomicSpeciesCard,AtomicPositionsCard,AtomicForcesCard}) =
    EachAtom(card)

Base.length(iter::EachAtom) = length(iter.card.data)

Base.iterate(iter::EachAtom, state=1) =
    state > length(iter) ? nothing : (iter.card.data[state], state + 1)

Base.eltype(iter::EachAtom) = eltype(iter.card.data)

"""
    optconvert(new_option::AbstractString, card::AbstractCellParametersCard)

Convert the option of an `AbstractCellParametersCard` from "bohr" to "angstrom", or its reverse.

!!! warning
    It does not support conversion between `"alat"` and the others.
"""
function optconvert(new_option::AbstractString, card::AbstractCellParametersCard)
    old_option = optionof(card)
    if new_option == old_option
        return card  # No conversion is needed
    else
        typeof(card)(
            if (old_option => new_option) == ("bohr" => "angstrom")
                return @. ustrip(u"angstrom", card.data * u"bohr")
            elseif (old_option => new_option) == ("angstrom" => "bohr")
                return @. ustrip(u"bohr", card.data * u"angstrom")
            else
                error("unknown conversion rule $(old_option => new_option)!")
            end,
        )
    end
end # function optconvert

abstract type KPointsCard <: Card end

struct KMeshCard <: KPointsCard
    data::MonkhorstPackGrid
end

struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard(data, option)

Represent the `K_POINTS` card in QE.

# Arguments
- `data::Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}`: A Γ point, a Monkhorst--Pack grid or a vector containing `SpecialKPoint`s.
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
"""
@struct_hash_equal struct SpecialPointsCard <: KPointsCard
    data::Vector{ReciprocalPoint}
    option::String
    function SpecialPointsCard(data, option="tpiba")
        @assert option in optionpool(SpecialPointsCard)
        return new(data, option)
    end
end
function SpecialPointsCard(data::AbstractMatrix, option="tpiba")
    @assert size(data, 2) == 4
    return SpecialPointsCard(map(x -> ReciprocalPoint(x...), eachrow(data)), option)
end

optionof(::KMeshCard) = "automatic"
optionof(::GammaPointCard) = "gamma"

optionpool(::Type{AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
optionpool(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
optionpool(::Type{KMeshCard}) = ("automatic",)
optionpool(::Type{GammaPointCard}) = ("gamma",)
optionpool(::Type{SpecialPointsCard}) =
    ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

groupname(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
groupname(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
groupname(::Type{CellParametersCard}) = "CELL_PARAMETERS"
groupname(::Type{<:KPointsCard}) = "K_POINTS"
