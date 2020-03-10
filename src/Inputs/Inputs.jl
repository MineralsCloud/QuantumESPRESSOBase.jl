"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using LinearAlgebra: det

using Compat: isnothing
using Crystallography: BravaisLattice
using Kaleido: @batchlens
using OrderedCollections: OrderedDict

using QuantumESPRESSOBase: InputEntry, titleof, to_qe
using QuantumESPRESSOBase.Setters: CellParametersSetter, LensMaker

import Crystallography
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Setters

export Card
export to_dict, dropdefault, getnamelists, getcards, compulsory_namelists, compulsory_cards, getoption, allowed_options, entryname, titleof

"""
    InputEntry

Represent any component of a `QuantumESPRESSOInput`.

Hierachy of `InputEntry`:
```
QuantumESPRESSOBase.InputEntry
├─ QuantumESPRESSOBase.Cards.Card
│  ├─ CP.AbstractCellParametersCard
│  │  ├─ CP.CellParametersCard
│  │  └─ CP.RefCellParametersCard
│  ├─ CP.AtomicForcesCard
│  ├─ CP.AtomicPositionsCard
│  ├─ CP.AtomicSpeciesCard
│  ├─ CP.AtomicVelocitiesCard
│  ├─ CP.KPointsCard
│  ├─ PHonon.AbstractCellParametersCard
│  │  └─ PHonon.CellParametersCard
│  ├─ PHonon.AtomicForcesCard
│  ├─ PHonon.AtomicPositionsCard
│  ├─ PHonon.AtomicSpeciesCard
│  ├─ PHonon.KPointsCard
│  ├─ PWscf.AbstractCellParametersCard
│  │  └─ PWscf.CellParametersCard
│  ├─ PWscf.AtomicForcesCard
│  ├─ PWscf.AtomicPositionsCard
│  ├─ PWscf.AtomicSpeciesCard
│  └─ PWscf.KPointsCard
└─ QuantumESPRESSOBase.Namelists.Namelist
   ├─ CP.CellNamelist
   ├─ CP.ControlNamelist
   ├─ CP.ElectronsNamelist
   ├─ CP.IonsNamelist
   ├─ CP.PressAiNamelist
   ├─ CP.SystemNamelist
   ├─ CP.WannierNamelist
   ├─ PHonon.DynmatNamelist
   ├─ PHonon.MatdynNamelist
   ├─ PHonon.PhNamelist
   ├─ PHonon.Q2rNamelist
   ├─ PWscf.BandsNamelist
   ├─ PWscf.CellNamelist
   ├─ PWscf.ControlNamelist
   ├─ PWscf.DosNamelist
   ├─ PWscf.ElectronsNamelist
   ├─ PWscf.IonsNamelist
   └─ PWscf.SystemNamelist
```
"""
abstract type InputEntry end

"""
    entryname(::Type{<:InputEntry})
    entryname(::InputEntry)

Return the field name of a `Namelist` or a `Card` in a `QuantumESPRESSOInput`.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> entryname(SystemNamelist) == entryname(SystemNamelist()) == :system
true
```
"""
entryname(x::InputEntry) = entryname(typeof(x))

"""
    titleof(::Type{<:InputEntry})
    titleof(::InputEntry)

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> titleof(SystemNamelist()) == titleof(SystemNamelist) == "SYSTEM"
true
```
"""
titleof(x::InputEntry) = titleof(typeof(x))

"""
    Namelist <: InputEntry

The abstraction of an component of a `QuantumESPRESSOInput`, a basic Fortran data structure.
"""
abstract type Namelist <: InputEntry end

"""
    to_dict(nml; defaultorder = true)

Convert a `Namelist` to a dictionary.

# Arguments
- `nml::Namelist`: the namelist to be converted.
- `defaultorder::Bool = true`: whether or not use the default order of parameters in QE's docs.
"""
function to_dict(nml::Namelist; defaultorder::Bool = true)
    dict = (defaultorder ? OrderedDict{Symbol,Any}() : Dict{Symbol,Any}())
    for n in propertynames(nml)
        dict[n] = getproperty(nml, n)
    end
    return dict
end

"""
    dropdefault(nml::Namelist)

Return an `AbstractDict` of non-default values of a `Namelist`.
"""
function dropdefault(nml::Namelist)
    default = typeof(nml)()  # Create a `Namelist` with all default values
    # Compare `default` with `nml`, discard the same values
    result = filter!(item -> item.second != getfield(default, item.first), to_dict(nml))
    isempty(result) && @info "Every entry in the namelist is the default value!"
    return result
end

"""
    Card <: InputEntry

The abstraction of all components of a `QuantumESPRESSOInput` that is not a `Namelist`.
"""
abstract type Card <: InputEntry end

"""
    getoption(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `getoption(x)` is suggested.
"""
getoption(card::Card) = getfield(card, :option)

"""
    allowed_options(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards, QuantumESPRESSOBase.Cards.PWscf

julia> allowed_options(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> allowed_options(CellParametersCard)
("alat", "bohr", "angstrom")

julia> allowed_options(KPointsCard)
("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
allowed_options(::Type{<:Card}) = nothing

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(f, obj) = Iterators.filter(f, (getfield(obj, i) for i in 1:nfields(obj)))

"""
    getnamelists(input::QuantumESPRESSOInput)

Return an iterable of `Namelist`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
getnamelists(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Namelist), input)

"""
    getcards(input::QuantumESPRESSOInput)

Return an iterable of `Card`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
getcards(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Card), input)

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

using .PWscf: PWInput
using .CP: CPInput

"""
    compulsory_namelists(input::Union{PWInput,CPInput})

Return an iterable of compulsory `Namelist`s of a `PWInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
It is lazy, you may want to `collect` it.

To get the optional `Namelist`s, use `(!compulsory_namelists)(input)` (Note the parenthesis!).
"""
compulsory_namelists(input::Union{PWInput,CPInput}) =
    (getfield(input, x) for x in (:control, :system, :electrons))
Base.:!(::typeof(compulsory_namelists)) =
    function (input::T) where {T<:Union{PWInput,CPInput}}
        (
            getfield(input, x) for x in fieldnames(T) if x ∉ (
                :control,
                :system,
                :electrons,
            ) && fieldtype(T, x) <: Namelist
        )
    end

"""
    compulsory_cards(input::PWInput)

Return an iterable of compulsory `Card`s of a `PWInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)` (Note the parenthesis!).
"""
compulsory_cards(input::PWInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions, :k_points))
"""
    compulsory_cards(input::CPInput)

Return an iterable of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)` (Note the parenthesis!).
"""
compulsory_cards(input::CPInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions))
Base.:!(::typeof(compulsory_cards)) = function (input::T) where {T<:Union{PWInput,CPInput}}
    (
        getfield(input, x) for x in fieldnames(T) if x ∉ (
            :atomic_species,
            :atomic_positions,
        ) && fieldtype(T, x) <: Card
    )
end

"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function Crystallography.cellvolume(input::PWInput)
    if isnothing(input.cell_parameters)
        return abs(det(BravaisLattice(input.system)()))
    else
        if getoption(input.cell_parameters) == "alat"
            # If no value of `celldm` is changed...
            isnothing(input.system.celldm[1]) && error("`celldm[1]` is not defined!")
            return input.system.celldm[1]^3 * abs(det(input.cell_parameters.data))
        else  # "bohr" or "angstrom"
            return cellvolume(input.cell_parameters)
        end
    end
end # function Crystallography.cellvolume

function QuantumESPRESSOBase.to_qe(
    nml::Namelist;
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
)
    namelist_name = titleof(nml)
    content = to_qe(dropdefault(nml); indent = indent, delim = delim)
    return "&$namelist_name" * newline * content * '/'
end
function QuantumESPRESSOBase.to_qe(
    input::QuantumESPRESSOInput;
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
    verbose::Bool = false,
)::String
    content = ""
    for namelist in getnamelists(input)
        content *=
            to_qe(namelist, indent = indent, delim = delim, newline = newline) * newline
    end
    for card in getcards(input)
        content *= to_qe(card, indent = indent, delim = delim, newline = newline) * newline
    end
    return content
end

function Setters.make(::LensMaker{CellParametersSetter,<:Union{PWInput,CPInput}})
    return @batchlens begin
        _.cell_parameters
        _.system.ibrav
        _.system.celldm
    end
end # function Setters.make

function Setters.preset_values(::CellParametersSetter, template::Union{PWInput,CPInput})
    # !isnothing(template.cell_parameters) && return template
    system = template.system
    return (
        CellParametersCard("alat", BravaisLattice(system)()),
        0,
        [system.celldm[1]],
    )
end # function Setters.preset_values

end
