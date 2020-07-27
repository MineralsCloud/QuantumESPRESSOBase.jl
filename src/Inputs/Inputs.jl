"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using AbInitioSoftwareBase.Inputs: Input
using Compat: only, isnothing
using Crystallography: Bravais, CellParameters, PrimitiveTriclinic
using OptionalArgChecks: @argcheck
using PyFortran90Namelists: fstring

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import AbInitioSoftwareBase.Inputs.Formats: delimiter, newline, indent
import Crystallography

export getoption,
    optionpool,
    titleof,
    inputstring,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    allnamelists,
    allcards

abstract type QuantumESPRESSOInputEntry end

"""
    Namelist <: QuantumESPRESSOInputEntry

The abstraction of an component of a `Input`, a basic Fortran data structure.
"""
abstract type Namelist <: QuantumESPRESSOInputEntry end

"""
    Card <: QuantumESPRESSOInputEntry

The abstraction of all components of a `Input` that is not a `Namelist`.
"""
abstract type Card <: QuantumESPRESSOInputEntry end

"""
    titleof(::Union{Namelist,Card})

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
titleof(x::QuantumESPRESSOInputEntry) = titleof(typeof(x))

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
    optionpool(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards, QuantumESPRESSOBase.Cards.PWscf

julia> optionpool(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> optionpool(CellParametersCard)
("alat", "bohr", "angstrom")

julia> optionpool(KPointsCard)
("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
function optionpool end

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput <: Input end

# This is a helper function and should not be exported.
entryname(S::Type{<:QuantumESPRESSOInputEntry}, T::Type{<:QuantumESPRESSOInput}) =
    only(fieldname(T, i) for (i, m) in enumerate(fieldtypes(T)) if S <: m)

function allnamelists end

function allcards end

function required_namelists end

function optional_namelists end

function required_cards end

function optional_cards end

# Do not export this type!
struct _Celldm{T<:Bravais}
    data::Any
    function _Celldm{T}(data) where {T}
        @argcheck 1 <= length(data) <= 6
        return new(data)
    end
end

function Base.getindex(x::_Celldm, i::Integer)
    a = x.data[1]
    if i == 1
        return a
    elseif i in 2:3
        return a * x.data[i]
    elseif i in 4:6
        return acos(x.data[10-i])
    else
        throw(BoundsError(x.data, i))
    end
end # function Base.getindex
function Base.getindex(x::_Celldm{PrimitiveTriclinic}, i::Integer)
    a = x.data[1]
    if i == 1
        return a
    elseif i in 2:3
        return a * x.data[i]
    elseif i in 4:6
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
# include("CP/CP.jl")
include("PHonon.jl")

"""
    inputstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function inputstring(input::QuantumESPRESSOInput)
    return join(
        map(Iterators.filter(!isnothing, getfield(input, i) for i = 1:nfields(input))) do f
            inputstring(f)
        end,
        newline(input),
    )
end
"""
    inputstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function inputstring(nml::Namelist)
    content = _nmlinputstring(
        dropdefault(nml);
        indent = indent(nml),
        delimiter = delimiter(nml),
        newline = newline(nml),
    )
    return join(filter(!isempty, ("&" * titleof(nml), content, '/')), newline(nml))
end
inputstring(x::AbstractString) = string(x)
function _nmlinputstring(
    dict::AbstractDict;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        map(keys(dict), values(dict)) do key, value
            _nmlinputstring(
                key,
                value;
                indent = indent,
                delimiter = delimiter,
                newline = newline,
            )
        end,
        newline,
    )
end
function _nmlinputstring(
    key,
    value::AbstractVector;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        map(Iterators.filter(x -> x[2] !== nothing, enumerate(value))) do (i, x)
            indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter)
        end,
        newline,
    )
end
function _nmlinputstring(
    key,
    value::NamedTuple;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        map(keys(value), values(value)) do x, y
            indent * join((string(key, '%', x), "=", fstring(y)), delimiter)
        end,
        newline,
    )
end
_nmlinputstring(key, value; indent = ' '^4, delimiter = ' ', newline = '\n') =
    indent * join((string(key), "=", fstring(value)), delimiter)

newline(::Union{QuantumESPRESSOInput,Namelist,Card}) = '\n'

indent(::Namelist) = ' '^4

delimiter(::Namelist) = ' '

end
