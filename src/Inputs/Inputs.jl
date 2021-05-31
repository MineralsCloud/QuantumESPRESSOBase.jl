"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using AbInitioSoftwareBase.Inputs: Input, InputEntry, Namelist, Card, Setter
using Compat: only, isnothing
using OptionalArgChecks: @argcheck
using PyFortran90Namelists: fstring

import AbInitioSoftwareBase.Inputs: asstring, groupname

export optionof,
    optionpool,
    groupname,
    asstring,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    allnamelists,
    allcards


"""
    groupname(::Union{Namelist,Card})

Return the title of the input entry in Quantum ESPRESSO.

The definition `groupname(x) = groupname(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Inputs.PWscf: ControlNamelist

julia> groupname(ControlNamelist()) == groupname(ControlNamelist) == "CONTROL"
true
```
"""
groupname(x::InputEntry) = groupname(typeof(x))

"""
    dropdefault(nml::Namelist)

Return an `AbstractDict` of non-default values of a `Namelist`.
"""
function dropdefault(nml::Namelist)
    default = typeof(nml)()  # Create a `Namelist` with all default values
    # Compare `default` with `nml`, discard the same values
    result = filter!(item -> item.second != getfield(default, item.first), Dict(nml))
    # for (drivingarg, passivearg) in _coupledargs(typeof(nml))
    # rule
    # end
    if isempty(result)
        @info "Every entry in the namelist is the default value!"
    end
    return result
end

# _coupledargs(::Type{<:Namelist}) = ()

Base.Dict(nml::Namelist) =
    Dict(name => getproperty(nml, name) for name in propertynames(nml))
Base.NamedTuple(nml::Namelist) =
    NamedTuple{propertynames(nml)}(getproperty(nml, name) for name in propertynames(nml))
Base.setdiff(a::T, b::T) where {T<:Namelist} = setdiff(Dict(a), Dict(b))

"""
    optionof(x::Card)

Return a `String` representing the option for `Card` `x`.

!!! warning
    Do not use `x.option` to access a `Card`'s `option`.
"""
optionof(card::Card) = hasfield(typeof(card), :option) ? getfield(card, :option) : nothing

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

julia> optionpool(SpecialKPointsCard)
("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
function optionpool end

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput <: Input end

# This is a helper function and should not be exported.
entryname(S::Type{<:InputEntry}, T::Type{<:QuantumESPRESSOInput}) =
    only(fieldname(T, i) for (i, m) in enumerate(fieldtypes(T)) if S <: m)

function allnamelists end

function allcards end

function required_namelists end

function optional_namelists end

function required_cards end

function optional_cards end

struct VerbositySetter <: Setter
    v::String
    function VerbositySetter(v)
        @assert v in ("high", "low")
        return new(v)
    end
end

include("PWscf/PWscf.jl")
# include("CP/CP.jl")
include("PHonon/PHonon.jl")

"""
    asstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function asstring(input::QuantumESPRESSOInput)
    return join(
        map(
            asstring,
            Iterators.filter(!isnothing, getfield(input, i) for i in 1:nfields(input)),
        ),
        newline(input),
    ) * newline(input)  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    content = _nmlasstring(
        dropdefault(nml);
        indent = indent(nml),
        delimiter = delimiter(nml),
        newline = newline(nml),
    )
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), newline(nml))
end
asstring(x::AbstractString) = string(x)
function _nmlasstring(
    dict::AbstractDict;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        map(keys(dict), values(dict)) do key, value
            _nmlasstring(
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
function _nmlasstring(
    key,
    value::AbstractVector;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        (
            indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter) for
            (i, x) in enumerate(value) if !isnothing(x)
        ),
        newline,
    )
end
function _nmlasstring(
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
_nmlasstring(key, value; indent = ' '^4, delimiter = ' ', newline = '\n') =
    indent * join((string(key), "=", fstring(value)), delimiter)

newline(::Union{QuantumESPRESSOInput,Namelist,Card}) = '\n'

indent(::Namelist) = ' '^4

delimiter(::Namelist) = ' '

include("crystallography.jl")

end
