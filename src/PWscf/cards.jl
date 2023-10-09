using ConstructionBase: constructorof
using CrystallographyBase: Cell, MonkhorstPackGrid
using StaticArrays: MVector, MMatrix

import AbInitioSoftwareBase: listpotentials
import ..QuantumESPRESSOBase: SpecialPoint, getoption, optionpool

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    AtomicForce,
    AtomicForcesCard,
    AtomicVelocity,
    AtomicVelocitiesCard,
    KPointsCard,
    KMeshCard,
    GammaPointCard,
    SpecialPointsCard,
    CellParametersCard
export getoption, convertoption, optionpool, eachatom, listpotentials

"""
    AtomicSpecies(atom::Union{Char,String}, mass, pseudopot)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)

Represent each line in the `ATOMIC_SPECIES` card in `pw.x` input files.

The `atom` field accepts no more than 3 characters.

# Examples
```jldoctest
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
    function AtomicSpecies(atom, mass, pseudopot)
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
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
    listpotentials(card::AtomicSpeciesCard)

List the pseudopotentials in an `AtomicSpeciesCard`.
"""
listpotentials(card::AtomicSpeciesCard) = map(atom -> atom.pseudopot, eachatom(card))

"""
    AtomicPosition(atom::Union{Char,String}, pos[, if_pos])
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Represent each line in the `ATOMIC_POSITIONS` card in `pw.x` input files.

The `atom` field accepts no more than 3 characters.

# Examples
```jldoctest
julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])

julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```
"""
@struct_hash_equal struct AtomicPosition
    "Label of the atom as specified in `AtomicSpecies`."
    atom::String
    "Atomic positions. A three-element vector of floats."
    pos::MVector{3,Float64}
    """
    Component `i` of the force for this atom is multiplied by `if_pos(i)`,
    which must be either `0` or `1`.  Used to keep selected atoms and/or
    selected components fixed in MD dynamics or structural optimization run.

    With `crystal_sg` atomic coordinates the constraints are copied in all equivalent
    atoms.
    """
    if_pos::MVector{3,Bool}
    function AtomicPosition(atom, pos, if_pos=trues(3))
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
        return new(string(atom), pos, if_pos)
    end
end
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicPosition, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)

"""
    AtomicPositionsCard <: Card

Represent the `ATOMIC_POSITIONS` card in `pw.x` input files.

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
AtomicPositionsCard(cell::Cell, option="alat") = AtomicPositionsCard(
    map(cell.atoms, cell.positions) do atom, position
        AtomicPosition(string(atom), position)
    end,
    option,
)

"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end
"""
    CellParametersCard(data::AbstractMatrix, option::String)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
struct CellParametersCard <: AbstractCellParametersCard
    data::MMatrix{3,3,Float64,9}
    option::String
    function CellParametersCard(data, option="alat")
        @assert option in optionpool(CellParametersCard)
        return new(data, option)
    end
end
CellParametersCard(lattice::Lattice, option="alat") =
    CellParametersCard(transpose(lattice.data), option)
CellParametersCard(lattice::Lattice{<:Length}) =
    CellParametersCard(Lattice(map(x -> ustrip(u"bohr", x), lattice.data)), "bohr")
CellParametersCard(cell::Cell, option="alat") =
    CellParametersCard(transpose(cell.lattice), option)

@struct_hash_equal struct AtomicForce
    atom::String
    force::MVector{3,Float64}
    function AtomicForce(atom, force)
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
        return new(string(atom), force)
    end
end

@struct_hash_equal struct AtomicForcesCard <: Card
    data::Vector{AtomicForce}
end

@struct_hash_equal struct AtomicVelocity
    atom::String
    velocity::MVector{3,Float64}
    function AtomicVelocity(atom, velocity)
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
        return new(string(atom), velocity)
    end
end

@struct_hash_equal struct AtomicVelocitiesCard <: Card
    data::Vector{AtomicVelocity}
end

# See https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1008-L1032 & https://github.com/JuliaLang/julia/blob/de3a70a/base/io.jl#L971-L1054
struct EachAtom{T}
    card::T
end

eachatom(card::Card) = EachAtom(card)

Base.length(iter::EachAtom) = length(iter.card.data)

Base.iterate(iter::EachAtom, state=1) = iterate(iter.card.data, state)

Base.eltype(iter::EachAtom) = eltype(iter.card.data)

"""
    convertoption(card::AbstractCellParametersCard, new_option::AbstractString)

Convert the option of an `AbstractCellParametersCard` from "bohr" to "angstrom", etc.

!!! warning
    It does not support conversions between `"alat"` and others.
"""
function convertoption(card::AbstractCellParametersCard, new_option::AbstractString)
    old_option = getoption(card)
    if new_option == old_option
        return card  # No conversion is needed
    else
        constructor = constructorof(typeof(card))
        if (old_option => new_option) == ("bohr" => "angstrom")
            return constructor(ustrip.(u"angstrom", card.data .* u"bohr"))
        elseif (old_option => new_option) == ("angstrom" => "bohr")
            return constructor(ustrip.(u"bohr", card.data .* u"angstrom"))
        else
            error("unknown conversion rule $(old_option => new_option)!")
        end
    end
end

"Represent the general `K_POINTS` card in Quantum ESPRESSO."
abstract type KPointsCard <: Card end
"""
    KMeshCard(data::MonkhorstPackGrid)

Represent the `K_POINTS` card in Quantum ESPRESSO with Monkhorst–Pack grid.
"""
struct KMeshCard <: KPointsCard
    data::MonkhorstPackGrid
end
"""
    GammaPointCard()

Represent the `K_POINTS` card in Quantum ESPRESSO with only Γ-point.
"""
struct GammaPointCard <: KPointsCard end
"""
    SpecialPointsCard(data, option)

Represent the `K_POINTS` card in Quantum ESPRESSO with a list of k-points.

# Arguments
- `data::Vector{SpecialPoint}`: a vector containing `SpecialPoint`s.
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
"""
@struct_hash_equal struct SpecialPointsCard <: KPointsCard
    data::Vector{SpecialPoint}
    option::String
    function SpecialPointsCard(data, option="tpiba")
        @assert option in optionpool(SpecialPointsCard)
        return new(data, option)
    end
end

struct AdditionalKPointsCard <: Card end

struct ConstraintsCard <: Card end

struct OccupationsCard <: Card end

struct SolventsCard <: Card end

struct HubbardCard <: Card end

getoption(::KMeshCard) = "automatic"
getoption(::GammaPointCard) = "gamma"

optionpool(card::Card) = optionpool(typeof(card))
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
