"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using LinearAlgebra: det

using Compat: isnothing, only
using Crystallography: BravaisLattice
using Kaleido: @batchlens
using OrderedCollections: OrderedDict

using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Setters: CellParametersSetter, LensMaker

import Crystallography
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Setters

export Card
export to_dict, dropdefault, getnamelists, getcards, getoption, allowed_options, titleof

"""
    InputEntry

Represent any component of a `QuantumESPRESSOInput`.

Hierachy of `InputEntry`:
```

```
"""
abstract type InputEntry end

"""
    titleof(::Type{<:InputEntry})
    titleof(::InputEntry)

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Inputs.PWscf: ControlNamelist

julia> titleof(ControlNamelist()) == titleof(ControlNamelist) == "CONTROL"
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

!!! warning
    A user should not use `x.option` to access a `Card`'s `option`.
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

# This is a helper function and should not be exported.
entryname(S::Type{<:InputEntry}, T::Type{<:QuantumESPRESSOInput}) = only(fieldname(T, i) for (i, m) in enumerate(fieldtypes(T)) if S <: m)

"""
    getnamelists(input::QuantumESPRESSOInput, selector::Symbol = :all)

Return an iterable of `Namelist`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.

Return an iterable of compulsory `Namelist`s of a `PWInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
"""
getnamelists(input::T, selector::Symbol = :all) where {T<:QuantumESPRESSOInput} = (getfield(input, x) for x in _selectnamelists(T, Val(selector)))

"""
    getcards(input::T, selector::Symbol = :all)

Return an iterable of `Card`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.

Return an iterable of compulsory `Card`s of a `PWInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
Return an iterable of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
"""
getcards(input::T, selector::Symbol = :all) where {T<:QuantumESPRESSOInput} = (getfield(input, x) for x in _selectcards(T, Val(selector)))

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

using .PWscf: PWInput
using .CP: CPInput

_selectnamelists(T::Type{<:QuantumESPRESSOInput}, ::Val{:all}) = Tuple(entryname(x, T) for x in fieldtypes(T) if x <: Namelist)
_selectnamelists(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:compulsory}) = (:control, :system, :electrons)
_selectnamelists(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:optional}) = setdiff(_selectnamelists(T, Val(:all)), _selectnamelists(T, Val(:compulsory)))

_selectcards(T::Type{<:QuantumESPRESSOInput}, ::Val{:all}) = Tuple(entryname(nonnothingtype(x), T) for x in fieldtypes(T) if x <: Union{Card,Nothing})
_selectcards(T::Type{PWInput}, ::Val{:compulsory}) = (:atomic_species, :atomic_positions, :k_points)
_selectcards(T::Type{CPInput}, ::Val{:compulsory}) = (:atomic_species, :atomic_positions)
_selectcards(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:optional}) = setdiff(_selectcards(T, Val(:all)), _selectcards(T, Val(:compulsory)))

# Referenced from https://discourse.julialang.org/t/how-to-get-the-non-nothing-type-from-union-t-nothing/30523
nonnothingtype(::Type{T}) where {T} = Core.Compiler.typesubtract(T, Nothing)  # Should not be exported

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
