module CP

using Compat: eachrow
using Formatting: sprintf1
using LinearAlgebra: det
using Parameters: @with_kw
using Setfield: get, set, @lens, @set
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful
using UnitfulAtomic

using ..Inputs:
    Namelist,
    QuantumESPRESSOInput,
    Card,
    optionof,
    optionpool,
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import Crystallography: Bravais, Lattice
# import Pseudopotentials: pseudoformat
import ..Inputs:
    optionpool,
    allnamelists,
    allcards,
    optionof,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards
using ..Formats: delimiter, newline, indent, floatfmt, intfmt

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    PressAiNamelist,
    WannierNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    RefCellParametersCard,
    AtomicVelocity,
    AtomicVelocitiesCard,
    AtomicForce,
    AtomicForcesCard,
    CPInput,
    optconvert

include("namelists.jl")
include("cards.jl")

"""
    CPInput <: QuantumESPRESSOInput
    CPInput(control, system, electrons, ions, cell, press_ai, wannier, atomic_species, atomic_positions, atomic_velocities, cell_parameters, ref_cell_parameters, constraints, occupations, atomic_forces, plot_wannier, autopilot)

Construct a `CPInput` which represents the input of program `cp.x`.

# Arguments
- `control::ControlNamelist=ControlNamelist()`: the `CONTROL` namelist of the input. Optional.
- `system::SystemNamelist=SystemNamelist()`: the `SYSTEM` namelist of the input. Optional.
- `electrons::ElectronsNamelist=ElectronsNamelist()`: the `ELECTRONS` namelist of the input. Optional.
- `ions::IonsNamelist=IonsNamelist()`: the `IONS` namelist of the input. Optional.
- `cell::CellNamelist=CellNamelist()`: the `CELL` namelist of the input. Optional.
- `press_ai::PressAiNamelist=PressAiNamelist()`: the `PRESS_AI` namelist of the input. Optional.
- `wannier::WannierNamelist=WannierNamelist()`: the `WANNIER` namelist of the input. Optional.
- `atomic_species::AtomicSpeciesCard`: the `ATOMIC_SPECIES` card of the input. Must be provided explicitly.
- `atomic_positions::AtomicPositionsCard`: the `ATOMIC_POSITIONS` card of the input. Must be provided explicitly.
- `atomic_velocities::AtomicVelocitiesCard`: the `ATOMIC_VELOCITIES` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
- `ref_cell_parameters::Union{Nothing,RefCellParametersCard}`: the `REF_CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
@with_kw struct CPInput <: QuantumESPRESSOInput
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    press_ai::PressAiNamelist = PressAiNamelist()
    wannier::WannierNamelist = WannierNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    atomic_velocities::AtomicVelocitiesCard
    cell_parameters::Union{Nothing,CellParametersCard} = nothing
    ref_cell_parameters::Union{Nothing,RefCellParametersCard} = nothing
    constraints::Union{Nothing,Float64} = nothing
    occupations::Union{Nothing,Float64} = nothing
    atomic_forces::Union{Nothing,AtomicForcesCard} = nothing
    plot_wannier::Union{Nothing,Float64} = nothing
    autopilot::Union{Nothing,Float64} = nothing
    @assert(
        !(cell_parameters === nothing && system.ibrav == 0),
        "Cannot specify `ibrav = 0` with an empty `cell_parameters`!"
    )
end # struct CPInput

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
        return typeof(card)(
            if (old_option => new_option) == ("bohr" => "angstrom")
                @. ustrip(u"angstrom", card.data * u"bohr")
            elseif (old_option => new_option) == ("angstrom" => "bohr")
                @. ustrip(u"bohr", card.data * u"angstrom")
            else
                error("unknown conversion rule $(old_option => new_option)!")
            end,
        )
    end
end # function optconvert

titleof(::Type{ControlNamelist}) = "CONTROL"
titleof(::Type{SystemNamelist}) = "SYSTEM"
titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
titleof(::Type{IonsNamelist}) = "IONS"
titleof(::Type{CellNamelist}) = "CELL"
titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"

"""
    Bravais(nml::SystemNamelist)

Return a `Bravais` from a `SystemNamelist`.
"""
Bravais(nml::SystemNamelist) = Bravais(nml.ibrav)

"""
    Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
Lattice(nml::SystemNamelist) = Lattice(Bravais(nml), nml.celldm)

"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))
"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function cellvolume(card::AbstractCellParametersCard)
    option = optionof(card)
    if option == "bohr"
        abs(det(card.data))
    elseif option == "angstrom"
        ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
        error("information not enough! Parameter `celldm[1]` needed!")
    end
end # function cellvolume

function inputstring(data::AtomicSpecies)
    return join(
        (sprintf1("%3s", data.atom), sprintf1(floatfmt(data), data.mass), data.pseudopot),
        delimiter(data),
    )
end
function inputstring(card::AtomicSpeciesCard)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES" *
           newline(card) *
           join((indent(card) * inputstring(x) for x in unique(card.data)), newline(card))
end
function inputstring(data::AtomicPosition)
    f(x) = x ? "" : "0"
    return join(
        [
            sprintf1("%3s", data.atom)
            map(x -> sprintf1(floatfmt(data), x), data.pos)
            map(f, data.if_pos)
        ],
        delimiter(data),
    )
end
function inputstring(card::AtomicPositionsCard)
    return "ATOMIC_POSITIONS { $(optionof(card)) }" *
           newline(card) *
           join((indent(card) * inputstring(x) for x in card.data), newline(card))
end
function inputstring(card::CellParametersCard)
    it = (
        indent * join((sprintf1(floatfmt(card), x) for x in row), delimiter(card)) for
        row in eachrow(card.data)
    )
    return "CELL_PARAMETERS { $(optionof(card)) }" * newline(card) * join(it, newline)
end

optionof(::AtomicVelocitiesCard) = "a.u"

optionpool(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
optionpool(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
optionpool(::Type{<:AtomicVelocity}) = ("a.u",)
optionpool(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

allnamelists(::Type{CPInput}) =
    (:control, :system, :electrons, :ions, :cell, :press_ai, :wannier)
allnamelists(x::CPInput) = (getfield(x, f) for f in allnamelists(typeof(x)))

allcards(::Type{CPInput}) = (
    :atomic_species,
    :atomic_positions,
    :atomic_velocities,
    :cell_parameters,
    :ref_cell_parameters,
    :constraints,
    :occupations,
    :atomic_forces,
    :plot_wannier,
    :autopilot,
)
allcards(x::CPInput) = (getfield(x, f) for f in allcards(typeof(x)))

required_namelists(::Type{CPInput}) = (:control, :system, :electrons)
required_namelists(x::CPInput) = (getfield(x, f) for f in required_namelists(typeof(x)))

optional_namelists(::Type{CPInput}) = (:ions, :cell, :press_ai, :wannier)
optional_namelists(x::CPInput) = (getfield(x, f) for f in optional_namelists(typeof(x)))

required_cards(::Type{CPInput}) = (:atomic_species, :atomic_positions)
required_cards(x::CPInput) = (getfield(x, f) for f in required_cards(typeof(x)))

optional_cards(::Type{CPInput}) = (
    :atomic_velocities,
    :cell_parameters,
    :ref_cell_parameters,
    :constraints,
    :occupations,
    :atomic_forces,
    :plot_wannier,
    :autopilot,
)
optional_cards(x::CPInput) = (getfield(x, f) for f in optional_cards(typeof(x)))

end
