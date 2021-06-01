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
using ConstructionBase: setproperties
using Crystallography: Cell
using Formatting: sprintf1
using LinearAlgebra: det, norm
using OptionalArgChecks: @argcheck
using Setfield: @set!
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful:
    AbstractQuantity, NoUnits, Temperature, dimension, upreferred, unit, ustrip, @u_str
using UnitfulAtomic

using ..Inputs: QuantumESPRESSOInput, Card, VerbositySetter, entryname

import AbInitioSoftwareBase.Inputs: InputEntry, Namelist, Setter, asstring, groupname
import AbInitioSoftwareBase.Inputs.Formatter: delimiter, newline, indent, floatfmt, intfmt
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
    VerbositySetter,
    ElectronicTemperatureSetter,
    ElecTempSetter,
    VolumeSetter,
    PressureSetter,
    StructureSetter,
    optconvert,
    xmldir,
    wfcfiles,
    exitfile,
    mkexitfile,
    optionof,
    optionpool,
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    asstring

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
    foreach(atomic_species.data) do datum
        path = joinpath(control.pseudo_dir, datum.pseudopot)
        if !isfile(path)
            @warn "pseudopotential file \"$path\" does not exist!"
        end
    end
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
PWInput(args::InputEntry...) = PWInput(; map(args) do arg
    entryname(typeof(arg), PWInput) => arg  # See https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/10
end...)

exitfile(template::PWInput) = abspath(
    expanduser(joinpath(template.control.outdir, template.control.prefix * ".EXIT")),
)
function mkexitfile(template::PWInput)
    path = exitfile(template)
    mkpath(dirname(path))
    return touch(path)
end

optionpool(::Type{AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
optionpool(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
optionpool(::Type{KMeshCard}) = ("automatic",)
optionpool(::Type{GammaPointCard}) = ("gamma",)
optionpool(::Type{SpecialPointsCard}) =
    ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

groupname(::Type{ControlNamelist}) = "CONTROL"
groupname(::Type{SystemNamelist}) = "SYSTEM"
groupname(::Type{ElectronsNamelist}) = "ELECTRONS"
groupname(::Type{IonsNamelist}) = "IONS"
groupname(::Type{CellNamelist}) = "CELL"
groupname(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
groupname(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
groupname(::Type{CellParametersCard}) = "CELL_PARAMETERS"
groupname(::Type{<:KPointsCard}) = "K_POINTS"

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

include("crystallography.jl")
include("asstring.jl")
include("set.jl")

end
