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
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful
using UnitfulAtomic

using ..Inputs:
    InputEntry,
    Namelist,
    QuantumESPRESSOInput,
    entryname,
    Card,
    _Celldm,
    getoption,
    allowed_options,
    inputstring

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import Crystallography
import Pseudopotentials
import ..Inputs

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
    wfcfiles

include("nml.jl")
include("card.jl")

Inputs.allowed_options(::Type{AtomicPositionsCard}) =
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

titleof(::Type{ControlNamelist}) = "CONTROL"
titleof(::Type{SystemNamelist}) = "SYSTEM"
titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
titleof(::Type{IonsNamelist}) = "IONS"
titleof(::Type{CellNamelist}) = "CELL"
titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"
titleof(::Type{<:KPointsCard}) = "K_POINTS"

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
    return "ATOMIC_SPECIES" *
           newline *
           join(map(unique(card.data)) do x
               indent * inputstring(x; delim = delim, numfmt = numfmt)
           end, newline)
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
           join(map(card.data) do x
               indent * inputstring(x; delim = delim, numfmt = numfmt)
           end, newline)
end
function inputstring(
    card::CellParametersCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    return "CELL_PARAMETERS { $(getoption(card)) }" *
           newline *
           join(map(eachrow(card.data)) do row
               indent * join((sprintf1(numfmt, x) for x in row), delim)
           end, newline)
end
inputstring(data::GammaPoint) = ""
function inputstring(data::MonkhorstPackGrid; delim = ' ', numfmt = "%5d", args...)
    return join(map([data.mesh; data.is_shift]) do x
        sprintf1(numfmt, x)
    end, delim)
end
inputstring(data::SpecialKPoint; delim = ' ', numfmt = "%14.9f", args...) =
    join(map(x -> sprintf1(numfmt, x), data), delim)
function inputstring(
    card::KPointsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    content = "K_POINTS { $(card.option) }" * newline
    if getoption(card) in ("gamma", "automatic")
        content *= indent * inputstring(card.data)
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= string(length(card.data), newline)
        content *= join(map(card.data) do x
            indent * inputstring(x; delim = delim, numfmt = numfmt)
        end, newline)
    end
    return content
end

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
    @assert(
        !(cell_parameters === nothing && system.ibrav == 0),
        "Cannot specify an empty `cell_parameters` with `ibrav = 0`!"
    )
end # struct PWInput
PWInput(args...) =
    PWInput(; Dict(zip(map(arg -> entryname(typeof(arg), PWInput), args), args))...)  # See https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/6

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
        abs(det(card.data))
    elseif option == "angstrom"
        ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
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

end
