using ConstructionBase: constructorof
using CrystallographyCore: EachAtom
using CrystallographyBase: Cell, MonkhorstPackGrid
using StaticArrays: FieldVector, SMatrix, Size

import AbInitioSoftwareBase: listpotentials
import CrystallographyCore: eachatom
import StaticArrays: similar_type

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
export getoption, convertoption, optionpool, listpotentials

abstract type AtomicData end

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
struct AtomicSpecies <: AtomicData
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

struct Position{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct IfPosition{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

# See https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.FieldVector
similar_type(::Type{<:Position}, ::Type{T}, s::Size{(3,)}) where {T} = Position{T}
similar_type(::Type{<:IfPosition}, ::Type{T}, s::Size{(3,)}) where {T} = IfPosition{T}

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
struct AtomicPosition <: AtomicData
    "Label of the atom as specified in `AtomicSpecies`."
    atom::String
    "Atomic positions. A three-element vector of floats."
    pos::Position{Float64}
    """
    Component `i` of the force for this atom is multiplied by `if_pos(i)`,
    which must be either `0` or `1`.  Used to keep selected atoms and/or
    selected components fixed in MD dynamics or structural optimization run.

    With `crystal_sg` atomic coordinates the constraints are copied in all equivalent
    atoms.
    """
    if_pos::IfPosition{Bool}
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
- `option::Symbol="alat"`: allowed values are: "alat", "bohr", "angstrom", "crystal", and "crystal_sg".
"""
@struct_hash_equal struct AtomicPositionsCard <: Card
    data::Vector{AtomicPosition}
    option::Symbol
    function AtomicPositionsCard(data, option=:alat)
        @assert option in optionpool(AtomicPositionsCard)
        return new(data, option)
    end
end
AtomicPositionsCard(cell::Cell, option=:alat) = AtomicPositionsCard(
    map(cell.atoms, cell.positions) do atom, position
        AtomicPosition(string(atom), position)
    end,
    option,
)

"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end
"""
    CellParametersCard(data::AbstractMatrix, option::Symbol)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
struct CellParametersCard <: AbstractCellParametersCard
    data::SMatrix{3,3,Float64,9}
    option::Symbol
    function CellParametersCard(data, option=:alat)
        @assert option in optionpool(CellParametersCard)
        return new(data, option)
    end
end
CellParametersCard(lattice::Lattice, option=:alat) =
    CellParametersCard(transpose(parent(lattice)), option)
CellParametersCard(lattice::Lattice{<:Length}) =
    CellParametersCard(Lattice(map(Base.Fix1(ustrip, u"bohr"), parent(lattice))), "bohr")
CellParametersCard(cell::Cell, option="alat") = CellParametersCard(Lattice(cell), option)

struct Force{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct Velocity{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

similar_type(::Type{<:Force}, ::Type{T}, s::Size{(3,)}) where {T} = Force{T}
similar_type(::Type{<:Velocity}, ::Type{T}, s::Size{(3,)}) where {T} = Velocity{T}

struct AtomicForce <: AtomicData
    atom::String
    force::Force{Float64}
    function AtomicForce(atom, force)
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
        return new(string(atom), force)
    end
end

@struct_hash_equal struct AtomicForcesCard <: Card
    data::Vector{AtomicForce}
end

struct AtomicVelocity <: AtomicData
    atom::String
    velocity::Velocity{Float64}
    function AtomicVelocity(atom, velocity)
        @assert length(atom) <= 3 "`atom` accepts no more than 3 characters!"
        return new(string(atom), velocity)
    end
end

@struct_hash_equal struct AtomicVelocitiesCard <: Card
    data::Vector{AtomicVelocity}
end

eachatom(card::AtomicPositionsCard) = EachAtom(
    Tuple(datum.atom for datum in card.data), Tuple(datum.pos for datum in card.data)
)

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
        if old_option == :bohr && new_option == :angstrom
            return constructor(ustrip.(u"angstrom", card.data .* u"bohr"))
        elseif old_option == :angstrom && new_option == :bohr
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
    option::Symbol
    function SpecialPointsCard(data, option=:tpiba)
        @assert option in optionpool(SpecialPointsCard)
        return new(data, option)
    end
end

struct AdditionalKPointsCard <: Card end

struct ConstraintsCard <: Card end

struct OccupationsCard <: Card end

struct SolventsCard <: Card end

struct HubbardCard <: Card end

getoption(::KMeshCard) = :automatic
getoption(::GammaPointCard) = :gamma

optionpool(card::Card) = optionpool(typeof(card))
optionpool(::Type{AtomicPositionsCard}) = (:alat, :bohr, :angstrom, :crystal, :crystal_sg)
optionpool(::Type{CellParametersCard}) = (:alat, :bohr, :angstrom)
optionpool(::Type{KMeshCard}) = (:automatic,)
optionpool(::Type{GammaPointCard}) = (:gamma,)
optionpool(::Type{SpecialPointsCard}) =
    (:tpiba, :crystal, :tpiba_b, :crystal_b, :tpiba_c, :crystal_c)

groupname(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
groupname(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
groupname(::Type{CellParametersCard}) = "CELL_PARAMETERS"
groupname(::Type{KMeshCard}) = "K_POINTS"
groupname(::Type{GammaPointCard}) = "K_POINTS"
groupname(::Type{SpecialPointsCard}) = "K_POINTS"
