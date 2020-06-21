module CP

using Compat: eachrow
using Crystallography: Bravais, Lattice, CellParameters, Cell
using Formatting: sprintf1
using LinearAlgebra: det
using Parameters: @with_kw
using Pseudopotentials: pseudopot_format
using Setfield: get, set, @lens, @set
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful
using UnitfulAtomic

using ..Inputs:
    Namelist,
    QuantumESPRESSOInput,
    Card,
    getoption,
    allowed_options,
    allnamelists,
    allcards,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import Crystallography
import Pseudopotentials
import ..Inputs:
    allowed_options,
    allnamelists,
    allcards,
    getoption,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards

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
    old_option = getoption(card)
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
    Crystallography.Bravais(nml::SystemNamelist)

Return a `Bravais` from a `SystemNamelist`.
"""
Crystallography.Bravais(nml::SystemNamelist) = Bravais(nml.ibrav)

"""
    Crystallography.Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
Crystallography.Lattice(nml::SystemNamelist) = Lattice(Bravais(nml), nml.celldm)

"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
Crystallography.cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))
"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function Crystallography.cellvolume(card::AbstractCellParametersCard)
    option = getoption(card)
    if option == "bohr"
        abs(det(card.data))
    elseif option == "angstrom"
        ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
        error("information not enough! Parameter `celldm[1]` needed!")
    end
end # function Crystallography.cellvolume

function inputstring(data::AtomicSpecies; delim = ' ', numfmt = "%14.9f", args...)
    return join(
        (sprintf1("%3s", data.atom), sprintf1(numfmt, data.mass), data.pseudopot),
        delim,
    )
end
function inputstring(
    card::AtomicSpeciesCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%20.10f",
    newline = '\n',
)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES" *
           newline *
           join(
               (
                   indent * inputstring(x; delim = delim, numfmt = numfmt)
                   for x in unique(card.data)
               ),
               newline,
           )
end
function inputstring(data::AtomicPosition; delim = ' ', numfmt = "%14.9f", args...)
    f(x) = x ? "" : "0"
    return join(
        [
            sprintf1("%3s", data.atom)
            map(x -> sprintf1(numfmt, x), data.pos)
            map(f, data.if_pos)
        ],
        delim,
    )
end
function inputstring(
    card::AtomicPositionsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    return "ATOMIC_POSITIONS { $(getoption(card)) }" *
           newline *
           join(
               (indent * inputstring(x; delim = delim, numfmt = numfmt) for x in card.data),
               newline,
           )
end
function inputstring(
    card::CellParametersCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    it = (
        indent * join((sprintf1(numfmt, x) for x in row), delim) for
        row in eachrow(card.data)
    )
    return "CELL_PARAMETERS { $(getoption(card)) }" * newline * join(it, newline)
end

getoption(::AtomicVelocitiesCard) = "a.u"

allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

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

compulsory_namelists(::Type{CPInput}) = (:control, :system, :electrons)
compulsory_namelists(x::CPInput) = (getfield(x, f) for f in compulsory_namelists(typeof(x)))

optional_namelists(::Type{CPInput}) = (:ions, :cell, :press_ai, :wannier)
optional_namelists(x::CPInput) = (getfield(x, f) for f in optional_namelists(typeof(x)))

compulsory_cards(::Type{CPInput}) = (:atomic_species, :atomic_positions)
compulsory_cards(x::CPInput) = (getfield(x, f) for f in compulsory_cards(typeof(x)))

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
