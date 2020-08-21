"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using AutoHashEquals: @auto_hash_equals
using Compat: eachrow, isnothing
using Crystallography: Cell
using Formatting: sprintf1
using LinearAlgebra: det, norm
using OptionalArgChecks: @argcheck
using Setfield: @set!
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful: AbstractQuantity, NoUnits, upreferred, unit, ustrip, @u_str
using UnitfulAtomic

using ..Inputs: QuantumESPRESSOInputEntry, Namelist, QuantumESPRESSOInput, entryname, Card

import AbInitioSoftwareBase.Inputs:
    inputstring, titleof, set_verbosity, set_elec_temp, set_press_vol, set_cell
import AbInitioSoftwareBase.Inputs.Formats: delimiter, newline, indent, floatfmt, intfmt
import Crystallography: Bravais, Lattice, cellvolume
# import Pseudopotentials: pseudoformat
import ..Inputs:
    optionpool,
    optionof,
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards
# _coupledargs

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
    SpecialPoint,
    KPointsCard,
    KMeshCard,
    GammaPointCard,
    SpecialPointsCard,
    PWInput,
    optconvert,
    xmldir,
    wfcfiles,
    optionof,
    optionpool,
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    set_verbosity,
    set_elec_temp,
    set_cell,
    set_press_vol,
    inputstring

include("namelists.jl")
include("cards.jl")

function iscompatible(system::SystemNamelist, cell_parameters::CellParametersCard)
    ibrav, celldm = system.ibrav, system.celldm
    if iszero(ibrav)
        if optionof(cell_parameters) in ("bohr", "angstrom")
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
- `k_points::AbstractKPointsCard`: the `K_POINTS` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
struct PWInput <: QuantumESPRESSOInput
    control::ControlNamelist
    system::SystemNamelist
    electrons::ElectronsNamelist
    ions::IonsNamelist
    cell::CellNamelist
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard}
    constraints::Union{Union{Nothing,Float64}}
    occupations::Union{Nothing,Float64}
    atomic_forces::Union{Nothing,AtomicForcesCard}
end # struct PWInput
function PWInput(;
    control = ControlNamelist(),
    system,
    electrons = ElectronsNamelist(),
    ions = IonsNamelist(),
    cell = CellNamelist(),
    atomic_species,
    atomic_positions,
    k_points,
    cell_parameters = nothing,
    constraints = nothing,
    occupations = nothing,
    atomic_forces = nothing,
)
    @argcheck !isnothing(cell_parameters) || system.ibrav != 0 "`cell_parameters` is empty with `ibrav = 0`!"
    return PWInput(
        control,
        system,
        electrons,
        ions,
        cell,
        atomic_species,
        atomic_positions,
        k_points,
        cell_parameters,
        constraints,
        occupations,
        atomic_forces,
    )
end
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
    set_elec_temp(system::PWInput, temperature::Union{Real,AbstractQuantity})

Return a modified `PWInput`, with finite temperature set.

!!! warning
    Can be used with(out) units. If no unit is given, "Ry" is chosen.
"""
function set_elec_temp(template::PWInput, temperature)
    @set! template.system = set_elec_temp(template.system, temperature)
    return template
end # function set_elec_temp

function set_press_vol(template::PWInput, pressure::Real, volume::Real)
    @set! template.cell.press = pressure
    factor = cbrt(volume / cellvolume(template))
    if isnothing(template.cell_parameters) || optionof(template.cell_parameters) == "alat"
        @set! template.system.celldm[1] *= factor
    else
        @set! template.system.celldm = zeros(6)
        @set! template.cell_parameters =
            optconvert("bohr", CellParametersCard(template.cell_parameters.data * factor))
    end
    return template
end # function set_press_vol
set_press_vol(template::PWInput, pressure::AbstractQuantity, volume::AbstractQuantity) =
    set_press_vol(template, ustrip(u"kbar", pressure), ustrip(u"bohr^3", volume))

function set_cell(template::PWInput, cell_parameters::CellParametersCard)
    if isnothing(template.cell_parameters)
        if optionof(cell_parameters) in ("bohr", "angstrom")
            @set! template.cell_parameters = cell_parameters
            @set! template.system.ibrav = 0
            @set! template.system.celldm = zeros(6)
        else
            @set! template.system.celldm = [template.system.celldm[1]]
            @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
        end
    else
        if optionof(template.cell_parameters) == "alat"
            if optionof(cell_parameters) in ("bohr", "angstrom")
                @set! template.system.celldm = [template.system.celldm[1]]
                cell_parameters = CellParametersCard(
                    cell_parameters.data / template.system.celldm[1],
                    "alat",
                )
            else  # "alat"
                @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
            end
        else
            if optionof(cell_parameters) == "alat"
                error("not matched!")
            end
        end
    end
    @set! template.cell_parameters = cell_parameters
    return template
end # function set_cell
function set_cell(template::PWInput, atomic_positions::AtomicPositionsCard)
    @set! template.atomic_positions = atomic_positions
    return template
end # function set_cell
set_cell(template::PWInput, c::CellParametersCard, a::AtomicPositionsCard) =
    set_cell(set_cell(template, c), a)
function set_cell(template::PWInput, cell::Cell, option1, option2)
    return set_cell(
        template,
        CellParametersCard(cell, option1),
        AtomicPositionsCard(cell, option2),
    )
end # function set_cell

optionpool(::Type{AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
optionpool(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
optionpool(::Type{KMeshCard}) = ("automatic",)
optionpool(::Type{GammaPointCard}) = ("gamma",)
optionpool(::Type{SpecialPointsCard}) =
    ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

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
    ("ATOMIC_POSITIONS { $(optionof(card)) }", map(inputstring, card.data)...),
    newline(card),
)
"""
    inputstring(card::CellParametersCard)

Return a `String` representing a `CellParametersCard`, valid for Quantum ESPRESSO's input.
"""
function inputstring(card::CellParametersCard)
    return join(
        (
            "CELL_PARAMETERS { $(optionof(card)) }",
            map(eachrow(card.data)) do row
                join((sprintf1(floatfmt(card), x) for x in row))
            end...,
        ),
        newline(card),
    )
end
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
inputstring(data::SpecialPoint) =
    indent(data) * join(map(x -> sprintf1(floatfmt(data), x), data), delimiter(data))
"""
    inputstring(card::KPointsCard)

Return a `String` representing a `KPointsCard`, valid for Quantum ESPRESSO's input.
"""
function inputstring(card::SpecialPointsCard)
    content = "K_POINTS { $(optionof(card)) }" * newline(card)
    return join((content, length(card.data), map(inputstring, card.data)...), newline(card))
end
inputstring(card::GammaPointCard) = "K_POINTS { $(optionof(card)) }" * newline(card)
function inputstring(card::KMeshCard)
    content = "K_POINTS { $(optionof(card)) }" * newline(card)
    return content * inputstring(card.data)
end

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
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function cellvolume(card::AbstractCellParametersCard)
    option = optionof(card)
    if option == "bohr"
        return abs(det(card.data))
    elseif option == "angstrom"
        return ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
        error("information not enough! Parameter `celldm[1]` needed!")
    end
end # function cellvolume
"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))
"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function cellvolume(input::PWInput)
    if isnothing(input.cell_parameters)
        return cellvolume(Lattice(input.system))
    else
        if optionof(input.cell_parameters) == "alat"
            # If no value of `celldm` is changed...
            if isnothing(input.system.celldm[1])
                error("`celldm[1]` is not defined!")
            else
                return input.system.celldm[1]^3 * abs(det(input.cell_parameters.data))
            end
        else  # "bohr" or "angstrom"
            return cellvolume(input.cell_parameters)
        end
    end
end # function cellvolume

"""
    allnamelists(x::PWInput)

Return an iterator of all `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
allnamelists(x::PWInput) = (getfield(x, f) for f in _allnamelists(typeof(x)))
_allnamelists(::Type{PWInput}) = (:control, :system, :electrons, :ions, :cell)

"""
    allcards(x::PWInput)

Get all `Card`s from a `PWInput`.
"""
allcards(x::PWInput) = (getfield(x, f) for f in _allcards(typeof(x)))
_allcards(::Type{PWInput}) = (
    :atomic_species,
    :atomic_positions,
    :k_points,
    :cell_parameters,
    :constraints,
    :occupations,
    :atomic_forces,
)

"""
    required_namelists(x::PWInput)

Return an iterator of required `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
required_namelists(x::PWInput) = (getfield(x, f) for f in _required_namelists(typeof(x)))
_required_namelists(::Type{PWInput}) = (:control, :system, :electrons)

"""
    optional_namelists(x::PWInput)

Return an iterator of optional `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
optional_namelists(x::PWInput) = (getfield(x, f) for f in _optional_namelists(typeof(x)))
_optional_namelists(::Type{PWInput}) = (:ions, :cell)

"""
    required_cards(x::PWInput)

Return an iterator of required `Card`s from a `PWInput`. You may want to `collect` them.
"""
required_cards(x::PWInput) = (getfield(x, f) for f in _required_cards(typeof(x)))
_required_cards(::Type{PWInput}) = (:atomic_species, :atomic_positions, :k_points)

"""
    optional_cards(x::PWInput)

Return an iterator of optional `Card`s from a `PWInput`. You may want to `collect` them.
"""
optional_cards(x::PWInput) = (getfield(x, f) for f in _optional_cards(typeof(x)))
_optional_cards(::Type{PWInput}) =
    (:cell_parameters, :constraints, :occupations, :atomic_forces)

indent(::Union{AtomicSpecies,AtomicPosition,SpecialPoint,MonkhorstPackGrid,AtomicForce}) =
    ' '^4

delimiter(
    ::Union{AtomicSpecies,AtomicPosition,SpecialPoint,MonkhorstPackGrid,AtomicForce},
) = ' '

floatfmt(::Union{AtomicSpecies,AtomicPosition,SpecialPoint}) = "%14.9f"
floatfmt(::CellParametersCard) = "%14.9f"

intfmt(::MonkhorstPackGrid) = "%5d"

end
