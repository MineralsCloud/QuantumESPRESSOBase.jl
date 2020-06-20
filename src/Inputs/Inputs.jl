"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using AbInitioSoftwareBase.Inputs: Input
using Compat: only
using Crystallography: Bravais, CellParameters, PrimitiveTriclinic
using Parameters: type2dict
using PyFortran90Namelists: fstring

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import Crystallography

export getoption,
    allowed_options,
    titleof,
    inputstring,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards,
    allnamelists,
    allcards

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

Base.Dict(nml::Namelist) =
    Dict(name => getproperty(nml, name) for name in propertynames(nml))
Base.NamedTuple(nml::Namelist) =
    (; (name => getproperty(nml, name) for name in propertynames(nml))...)
Base.setdiff(a::T, b::T) where {T<:Namelist} = setdiff(Dict(a), Dict(b))

"""
    getoption(x::Card)

Return the option for `Card` `x`.

!!! warning
    A user should not use `x.option` to access a `Card`'s `option`.
"""
getoption(card::Card) = hasfield(typeof(card), :option) ? getfield(card, :option) : nothing

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
function allowed_options end

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput <: Input end

# This is a helper function and should not be exported.
entryname(S::Type{<:InputEntry}, T::Type{<:QuantumESPRESSOInput}) =
    only(fieldname(T, i) for (i, m) in enumerate(fieldtypes(T)) if S <: m)

function allnamelists end

function allcards end

function compulsory_namelists end

function optional_namelists end

function compulsory_cards end

function optional_cards end

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

include("PWscf/PWscf.jl")
include("CP.jl")
include("PHonon.jl")

"""
    inputstring(input::QuantumESPRESSOInput; indent = ' '^4, delim = ' ', newline = '\n')

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function inputstring(input::QuantumESPRESSOInput; newline = '\n', kwargs...)
    return join(map(fieldnames(typeof(input))) do f
        x = getfield(input, f)
        if x !== nothing
            inputstring(x; newline = newline, kwargs...) * newline
        else
            ""
        end
    end)
end
"""
    inputstring(args::InputEntry...; indent = ' '^4, delim = ' ', newline = '\n')

Return a `String` representing a collection of `QuantumESPRESSOInput` fields, valid for Quantum ESPRESSO's input.
"""
function inputstring(args::InputEntry...; newline = '\n', kwargs...)
    return join(map(args) do x
        inputstring(x; newline = newline, kwargs...)
    end, newline)
end
"""
    inputstring(nml::Namelist; indent = ' '^4, delim = ' ', newline = '\n')

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function inputstring(nml::Namelist; newline = '\n', kwargs...)
    content = _inputstring(dropdefault(nml); newline = newline, kwargs...)
    return "&" * titleof(nml) * newline * content * '/'
end
function _inputstring(dict::AbstractDict; indent = ' '^4, delim = ' ', newline = '\n')
    return join(
        map(keys(dict), values(dict)) do key, value
            _inputstring(key, value; indent = indent, delim = delim, newline = newline)
        end,
    )
end
function _inputstring(
    key,
    value::AbstractVector;
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
)
    return join(
        map(enumerate(value)) do (i, x)
            if x !== nothing
                indent * join([string(key, '(', i, ')'), "=", fstring(x)], delim) * newline
            else
                ""
            end
        end,
    )
end
function _inputstring(key, value::NamedTuple; indent = ' '^4, delim = ' ', newline = '\n')
    return join(map(keys(value), values(value)) do x, y
        indent * join([string(key, '%', x), "=", fstring(y)], delim)
    end, newline)
end
_inputstring(key, value; indent = ' '^4, delim = ' ', newline = '\n') =
    indent * join([string(key), "=", fstring(value)], delim) * newline

end
