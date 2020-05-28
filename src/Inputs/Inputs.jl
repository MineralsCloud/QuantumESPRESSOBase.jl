"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using Compat: isnothing, only
using Crystallography: Bravais, Lattice, CellParameters, PrimitiveTriclinic, cellvolume
using Kaleido: @batchlens
using LinearAlgebra: det
using OrderedCollections: OrderedDict
using Parameters: type2dict
using PyFortran90Namelists: fstring

import Crystallography
import OrderedCollections
import ..Setters

export getnamelists, getcards, getoption, allowed_options, titleof, qestring

"""
    Namelist

The abstraction of an component of a `Input`, a basic Fortran data structure.
"""
abstract type Namelist end

"""
    Card

The abstraction of all components of a `Input` that is not a `Namelist`.
"""
abstract type Card end

const InputEntry = Union{Namelist,Card}

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
    dropdefault(nml::Namelist)

Return an `AbstractDict` of non-default values of a `Namelist`.
"""
function dropdefault(nml::Namelist)
    default = typeof(nml)()  # Create a `Namelist` with all default values
    # Compare `default` with `nml`, discard the same values
    result = filter!(item -> item.second != getfield(default, item.first), Dict(nml))
    isempty(result) && @info "Every entry in the namelist is the default value!"
    return result
end

Base.setdiff(a::T, b::T) where {T<:Namelist} = setdiff(type2dict(a), type2dict(b))

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
abstract type Input end

# This is a helper function and should not be exported.
entryname(S::Type{<:InputEntry}, T::Type{<:Input}) =
    only(fieldname(T, i) for (i, m) in enumerate(fieldtypes(T)) if S <: m)

"""
    getnamelists(input::Input, selector::Symbol = :all)

Return an iterable of `Namelist`s of a `Input`. It is lazy, you may want to `collect` it.

Return an iterable of compulsory `Namelist`s of a `PWInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
"""
getnamelists(input::T, selector::Symbol = :all) where {T<:Input} =
    (getfield(input, x) for x in _selectnamelists(T, Val(selector)))

"""
    getcards(input::T, selector::Symbol = :all)

Return an iterable of `Card`s of a `Input`. It is lazy, you may want to `collect` it.

Return an iterable of compulsory `Card`s of a `PWInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
Return an iterable of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
"""
getcards(input::T, selector::Symbol = :all) where {T<:Input} =
    (getfield(input, x) for x in _selectcards(T, Val(selector)))

# Do not export this type!
struct _Celldm{T<:Bravais}
    data::Any
    function _Celldm{T}(data) where {T}
        @assert 1 <= length(data) <= 6
        return new(data)
    end
end

function Base.getindex(x::_Celldm, i::Integer)
    a = x.data[1]
    if i == 1
        return a
    elseif i ∈ 2:3
        return a * x.data[i]
    elseif i ∈ 4:6
        return acos(x.data[10-i])
    else
        throw(BoundsError(x.data, i))
    end
end # function Base.getindex
function Base.getindex(x::_Celldm{PrimitiveTriclinic}, i::Integer)
    a = x.data[1]
    if i == 1
        return a
    elseif i ∈ 2:3
        return a * x.data[i]
    elseif i ∈ 4:6
        return acos(x.data[i])  # Note the difference!
    else
        throw(BoundsError(x.data, i))
    end
end # function Base.getindex
Base.getindex(x::_Celldm, I) = [x[i] for i in I]

function Base.convert(::Type{_Celldm{T}}, p::CellParameters) where {T}
    a, b, c, α, β, γ = p
    return _Celldm{T}([a, b / a, c / a, cos(γ), cos(β), cos(α)])  # What a horrible conversion!
end # function Base.convert
function Base.convert(::Type{_Celldm{PrimitiveTriclinic}}, p::CellParameters)
    a, b, c, α, β, γ = p
    return _Celldm{PrimitiveTriclinic}([a, b / a, c / a, cos(α), cos(β), cos(γ)])  # What a horrible conversion!
end # function Base.convert

Base.Dict(nml::Namelist) =
    Dict(name => getproperty(nml, name) for name in propertynames(nml))
OrderedCollections.OrderedDict(nml::Namelist) =
    OrderedDict(name => getproperty(nml, name) for name in propertynames(nml))

"""
    qestring(x; indent = ' '^4, delim = ' ')

Return a `String` representing the object, which is valid for Quantum ESPRESSO's input.
"""
function qestring(dict::AbstractDict; indent = ' '^4, delim = ' ')
    content = ""
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                isnothing(x) && continue
                content *= indent * join(["$key($i)", "=", "$(fstring(x))\n"], delim)
            end
        else
            content *= indent * join(["$key", "=", "$(fstring(value))\n"], delim)
        end
    end
    return content
end
qestring(::Nothing; args...) = ""
function qestring(nml::Namelist; indent = ' '^4, delim = ' ', newline = '\n')
    namelist_name = titleof(nml)
    content = qestring(dropdefault(nml); indent = indent, delim = delim)
    return "&$namelist_name" * newline * content * '/'
end
function qestring(input::Input; indent = ' '^4, delim = ' ', newline = '\n')
    content = ""
    for i in 1:nfields(input)
        content *=
            qestring(
                getfield(input, i),
                indent = indent,
                delim = delim,
                newline = newline,
            ) * newline
    end
    return content
end
function qestring(
    v::AbstractVector{<:InputEntry},
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
)
    content = ""
    for i in 1:length(v)
        content *=
            qestring(v[i], indent = indent, delim = delim, newline = newline) * newline
    end
    return content
end

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

using .PWscf: PWInput
using .CP: CPInput

_selectnamelists(T::Type{<:Input}, ::Val{:all}) =
    Tuple(entryname(x, T) for x in fieldtypes(T) if x <: Namelist)
_selectnamelists(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:compulsory}) =
    (:control, :system, :electrons)
_selectnamelists(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:optional}) =
    setdiff(_selectnamelists(T, Val(:all)), _selectnamelists(T, Val(:compulsory)))

_selectcards(T::Type{<:Input}, ::Val{:all}) = Tuple(
    entryname(nonnothingtype(x), T) for x in fieldtypes(T) if x <: Union{Card,Nothing}
)
_selectcards(T::Type{PWInput}, ::Val{:compulsory}) =
    (:atomic_species, :atomic_positions, :k_points)
_selectcards(T::Type{CPInput}, ::Val{:compulsory}) = (:atomic_species, :atomic_positions)
_selectcards(T::Union{Type{PWInput},Type{CPInput}}, ::Val{:optional}) =
    setdiff(_selectcards(T, Val(:all)), _selectcards(T, Val(:compulsory)))

# Referenced from https://discourse.julialang.org/t/how-to-get-the-non-nothing-type-from-union-t-nothing/30523
nonnothingtype(::Type{T}) where {T} = Core.Compiler.typesubtract(T, Nothing)  # Should not be exported

"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function Crystallography.cellvolume(input::PWInput)
    if isnothing(input.cell_parameters)
        return cellvolume(Lattice(input.system))
    else
        if getoption(input.cell_parameters) == "alat"
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
end # function Crystallography.cellvolume

end
