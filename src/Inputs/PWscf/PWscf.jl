"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: eachrow
using Crystallography: Bravais, Lattice, CellParameters, Cell, cellvolume
using Formatting: sprintf1
using LinearAlgebra: det
using Parameters: @with_kw
using Pseudopotentials: pseudopot_format
using Setfield: @set!
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful: AbstractQuantity, upreferred, unit, ustrip, @u_str
using UnitfulAtomic

using ..Inputs:
    QuantumESPRESSOInputEntry,
    Namelist,
    QuantumESPRESSOInput,
    entryname,
    Card,
    _Celldm,
    getoption

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import Crystallography
import Pseudopotentials
import ..Inputs:
    allowed_options,
    allnamelists,
    allcards,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards
import ..Formats: delimiter, newline, indent, floatfmt, intfmt

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    PWInput,
    optconvert,
    xmldir,
    wfcfiles,
    getoption,
    allowed_options,
    allnamelists,
    allcards,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards,
    set_verbosity,
    set_temperature,
    set_structure,
    inputstring

include("namelists.jl")
include("cards.jl")

function iscompatible(system::SystemNamelist, cell_parameters::CellParametersCard)
    ibrav, celldm = system.ibrav, system.celldm
    if iszero(ibrav)
        if getoption(cell_parameters) ∈ ("bohr", "angstrom")
            return all(iszero, celldm)
        else  # "alat"
            return !iszero(first(celldm))  # first(celldm) != 0
        end
    else
        return false
    end
end # function iscompatible
iscompatible(x::CellParametersCard, y::SystemNamelist) = iscompatible(y, x)

"""
    PWInput <: QuantumESPRESSOInput
    PWInput(control, system, electrons, ions, cell, atomic_species, atomic_positions, k_points, cell_parameters)

Construct a `PWInput` which represents the input of program `pw.x`.

# Arguments
- `control::ControlNamelist=ControlNamelist()`: the `CONTROL` namelist of the input. Optional.
- `system::SystemNamelist=SystemNamelist()`: the `SYSTEM` namelist of the input. Optional.
- `electrons::ElectronsNamelist=ElectronsNamelist()`: the `ELECTRONS` namelist of the input. Optional.
- `ions::IonsNamelist=IonsNamelist()`: the `IONS` namelist of the input. Optional.
- `cell::CellNamelist=CellNamelist()`: the `CELL` namelist of the input. Optional.
- `atomic_species::AtomicSpeciesCard`: the `ATOMIC_SPECIES` card of the input. Must be provided explicitly.
- `atomic_positions::AtomicPositionsCard`: the `ATOMIC_POSITIONS` card of the input. Must be provided explicitly.
- `k_points::KPointsCard`: the `K_POINTS` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
@with_kw struct PWInput <: QuantumESPRESSOInput
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard} = nothing
    constraints::Union{Union{Nothing,Float64}} = nothing
    occupations::Union{Nothing,Float64} = nothing
    atomic_forces::Union{Nothing,AtomicForcesCard} = nothing
    @assert cell_parameters !== nothing || system.ibrav != 0 "`cell_parameters` is empty with `ibrav = 0`!"
end # struct PWInput
PWInput(args::QuantumESPRESSOInputEntry...) = PWInput(; map(args) do arg
    entryname(typeof(arg), PWInput) => arg  # See https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/10
end...)

"""
    set_verbosity(template::PWInput, verbosity)

Return a modified `PWInput`, with verbosity set.
"""
function set_verbosity(template::PWInput, verbosity)
    @set! template.control = set_verbosity(template.control, verbosity)
    return template
end # function set_verbosity

"""
    set_temperature(system::PWInput, temperature)

Return a modified `PWInput`, with finite temperature set.

!!! warning
    Can be used with(out) units. If no unit is given, "Ry" is chosen.
"""
function set_temperature(template::PWInput, temperature)
    @set! template.system = set_temperature(template.system, temperature)
    return template
end # function set_temperature

function set_structure(
    template::PWInput,
    cell_parameters::Union{Nothing,CellParametersCard} = nothing,
    atomic_positions::Union{Nothing,AtomicPositionsCard} = nothing,
)
    if cell_parameters !== nothing
        @set! template.cell_parameters = cell_parameters
    end
    if atomic_positions !== nothing
        @set! template.atomic_positions = atomic_positions
    end
    return template
end # function set_structure
function set_structure(template::PWInput, cell::Cell, option1, option2)
    return set_structure(
        template,
        CellParametersCard(cell, option1),
        AtomicPositionsCard(cell, option2),
    )
end # function set_structure

allowed_options(::Type{AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{<:KPointsCard}) = (
    "tpiba",
    "automatic",
    "crystal",
    "gamma",
    "tpiba_b",
    "crystal_b",
    "tpiba_c",
    "crystal_c",
)

titleof(::Type{ControlNamelist}) = "CONTROL"
titleof(::Type{SystemNamelist}) = "SYSTEM"
titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
titleof(::Type{IonsNamelist}) = "IONS"
titleof(::Type{CellNamelist}) = "CELL"
titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
titleof(::Type{CellParametersCard}) = "CELL_PARAMETERS"
titleof(::Type{<:KPointsCard}) = "K_POINTS"

"""
    inputstring(data::AtomicSpecies)

Return a `String` representing a `AtomicSpecies`, valid for Quantum ESPRESSO's input.
"""
function inputstring(data::AtomicSpecies)
    return join(
        (
            indent(data),
            sprintf1("%3s", data.atom),
            sprintf1(floatfmt(data), data.mass),
            data.pseudopot,
        ),
        delimiter(data),
    )
end
"""
    inputstring(card::AtomicSpeciesCard)

Return a `String` representing a `AtomicSpeciesCard`, valid for Quantum ESPRESSO's input.
"""
inputstring(card::AtomicSpeciesCard) =
    join(("ATOMIC_SPECIES", map(inputstring, unique(card.data))...), newline(card))
"""
    inputstring(data::AtomicPosition)

Return a `String` representing a `AtomicPosition`, valid for Quantum ESPRESSO's input.
"""
function inputstring(data::AtomicPosition)
    content = join(
        (
            indent(data),
            sprintf1("%3s", data.atom),
            map(x -> sprintf1(floatfmt(data), x), data.pos)...,
        ),
        delimiter(data),
    )
    if !all(data.if_pos)
        return join((content, map(Int, data.if_pos)...), delimiter(data))
    else
        return content
    end
end
"""
    inputstring(card::AtomicPositionsCard)

Return a `String` representing a `AtomicPositionsCard`, valid for Quantum ESPRESSO's input.
"""
inputstring(card::AtomicPositionsCard) = join(
    ("ATOMIC_POSITIONS { $(getoption(card)) }", map(inputstring, card.data)...),
    newline(card),
)
"""
    inputstring(card::CellParametersCard)

Return a `String` representing a `CellParametersCard`, valid for Quantum ESPRESSO's input.
"""
function inputstring(card::CellParametersCard)
    return join(
        (
            "CELL_PARAMETERS { $(getoption(card)) }",
            map(eachrow(card.data)) do row
                join((sprintf1(floatfmt, x) for x in row))
            end...,
        ),
        newline(card),
    )
end
"""
    inputstring(data::GammaPoint)

Return a `String` representing a `GammaPoint`, valid for Quantum ESPRESSO's input.
"""
inputstring(data::GammaPoint) = indent(data)
"""
    inputstring(data::MonkhorstPackGrid)

Return a `String` representing a `MonkhorstPackGrid`, valid for Quantum ESPRESSO's input.
"""
function inputstring(data::MonkhorstPackGrid)
    return indent(data) * join(map([data.mesh; data.is_shift]) do x
        sprintf1(intfmt(data), x)
    end, delimiter(data))
end
"""
    inputstring(data::SpecialKPoint)

Return a `String` representing a `SpecialKPoint`, valid for Quantum ESPRESSO's input.
"""
inputstring(data::SpecialKPoint) =
    indent(data) * join(map(x -> sprintf1(floatfmt(data), x), data), delimiter(data))
"""
    inputstring(card::KPointsCard)

Return a `String` representing a `KPointsCard`, valid for Quantum ESPRESSO's input.
"""
function inputstring(card::KPointsCard)
    content = "K_POINTS { $(card.option) }"
    if getoption(card) in ("gamma", "automatic")
        return content * inputstring(card.data)
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        return join(
            (content, length(card.data), map(inputstring, card.data)...),
            newline(card),
        )
    end
end

"""
    Crystallography.Bravais(nml::SystemNamelist)

Return a `Bravais` from a `SystemNamelist`.
"""
Crystallography.Bravais(nml::SystemNamelist) = Bravais(nml.ibrav)

"""
    Crystallography.Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
function Crystallography.Lattice(nml::SystemNamelist)
    b = Bravais(nml)
    return Lattice(b, _Celldm{typeof(b)}(nml.celldm))
end # function Crystallography.Lattice

"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function Crystallography.cellvolume(card::AbstractCellParametersCard)
    option = getoption(card)
    if option == "bohr"
        return abs(det(card.data))
    elseif option == "angstrom"
        return ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
        error("information not enough! Parameter `celldm[1]` needed!")
    end
end # function Crystallography.cellvolume
"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
Crystallography.cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))
"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function Crystallography.cellvolume(input::PWInput)
    if input.cell_parameters === nothing
        return cellvolume(Lattice(input.system))
    else
        if getoption(input.cell_parameters) == "alat"
            # If no value of `celldm` is changed...
            if input.system.celldm[1] === nothing
                error("`celldm[1]` is not defined!")
            else
                return input.system.celldm[1]^3 * abs(det(input.cell_parameters.data))
            end
        else  # "bohr" or "angstrom"
            return cellvolume(input.cell_parameters)
        end
    end
end # function Crystallography.cellvolume

allnamelists(x::PWInput) = (getfield(x, f) for f in allnamelists(typeof(x)))
allnamelists(::Type{PWInput}) = (:control, :system, :electrons, :ions, :cell)

allcards(x::PWInput) = (getfield(x, f) for f in allcards(typeof(x)))
allcards(::Type{PWInput}) = (
    :atomic_species,
    :atomic_positions,
    :k_points,
    :cell_parameters,
    :constraints,
    :occupations,
    :atomic_forces,
)

compulsory_namelists(x::PWInput) = (getfield(x, f) for f in compulsory_namelists(typeof(x)))
compulsory_namelists(::Type{PWInput}) = (:control, :system, :electrons)

optional_namelists(x::PWInput) = (getfield(x, f) for f in optional_namelists(typeof(x)))
optional_namelists(::Type{PWInput}) = (:ions, :cell)

compulsory_cards(x::PWInput) = (getfield(x, f) for f in compulsory_cards(typeof(x)))
compulsory_cards(::Type{PWInput}) = (:atomic_species, :atomic_positions, :k_points)

optional_cards(x::PWInput) = (getfield(x, f) for f in optional_cards(typeof(x)))
optional_cards(::Type{PWInput}) =
    (:cell_parameters, :constraints, :occupations, :atomic_forces)

indent(
    ::Union{
        AtomicSpecies,
        AtomicPosition,
        GammaPoint,
        SpecialKPoint,
        MonkhorstPackGrid,
        AtomicForce,
    },
) = ' '^4

delimiter(
    ::Union{
        AtomicSpecies,
        AtomicPosition,
        GammaPoint,
        SpecialKPoint,
        MonkhorstPackGrid,
        AtomicForce,
    },
) = ' '

floatfmt(::Union{AtomicSpecies,AtomicPosition,SpecialKPoint}) = "%14.9f"
floatfmt(::CellParametersCard) = "%14.9f"

infmt(::MonkhorstPackGrid) = "%5d"

end
