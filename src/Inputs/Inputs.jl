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
using Crystallography:
    CellParameters,
    PrimitiveCubic,
    FaceCenteredCubic,
    BodyCenteredCubic,
    PrimitiveHexagonal,
    RCenteredHexagonal,
    PrimitiveTetragonal,
    BodyCenteredTetragonal,
    PrimitiveOrthorhombic,
    BCenteredOrthorhombic,
    ACenteredOrthorhombic,
    FaceCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    PrimitiveMonoclinic,
    CCenteredMonoclinic,
    BCenteredMonoclinic,
    PrimitiveTriclinic
using OptionalArgChecks: @argcheck
using PyFortran90Namelists: fstring

import AbInitioSoftwareBase.Inputs: inputstring, titleof
import AbInitioSoftwareBase.Inputs.Formats: delimiter, newline, indent
import Crystallography: Bravais, Lattice

export optionof,
    optionpool,
    titleof,
    inputstring,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    allnamelists,
    allcards

abstract type QuantumESPRESSOInputEntry end  # Define this to make the `eltype` not `Any` if both `Namelist` & `Card` exist.

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

Base.getindex(x::_Celldm, i) = getindex(x.data, i)
# function Base.getindex(x::_Celldm, i::Integer)
#     a = x.data[1]
#     if i == 1
#         return a
#     elseif i in 2:3
#         return a * x.data[i]
#     elseif i in 4:6
#         return acos(x.data[10-i])
#     else
#         throw(BoundsError(x.data, i))
#     end
# end # function Base.getindex
# function Base.getindex(x::_Celldm{PrimitiveTriclinic}, i::Integer)
#     a = x.data[1]
#     if i == 1
#         return a
#     elseif i in 2:3
#         return a * x.data[i]
#     elseif i in 4:6
#         return acos(x.data[i])  # Note the difference!
#     else
#         throw(BoundsError(x.data, i))
#     end
# end # function Base.getindex
# Base.getindex(x::_Celldm, I) = [x[i] for i in I]

# function Base.convert(::Type{_Celldm{T}}, p::CellParameters) where {T}
#     a, b, c, α, β, γ = p
#     return _Celldm{T}([a, b / a, c / a, cos(γ), cos(β), cos(α)])  # What a horrible conversion!
# end # function Base.convert
# function Base.convert(::Type{_Celldm{PrimitiveTriclinic}}, p::CellParameters)
#     a, b, c, α, β, γ = p
#     return _Celldm{PrimitiveTriclinic}([a, b / a, c / a, cos(α), cos(β), cos(γ)])  # What a horrible conversion!
# end # function Base.convert

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

function Bravais(ibrav::Integer)
    if ibrav == 1
        return PrimitiveCubic()
    elseif ibrav == 2
        return FaceCenteredCubic()
    elseif ibrav == 3
        return BodyCenteredCubic()
    elseif ibrav == 4
        return PrimitiveHexagonal()
    elseif ibrav == 5
        return RCenteredHexagonal()
    elseif ibrav == -5
        return RCenteredHexagonal()
    elseif ibrav == 6
        return PrimitiveTetragonal()
    elseif ibrav == 7
        return BodyCenteredTetragonal()
    elseif ibrav == 8
        return PrimitiveOrthorhombic()
    elseif ibrav == 9
        return BCenteredOrthorhombic()
    elseif ibrav == -9
        return BCenteredOrthorhombic()
    elseif ibrav == 91
        return ACenteredOrthorhombic()  # New in QE 6.5
    elseif ibrav == 10
        return FaceCenteredOrthorhombic()
    elseif ibrav == 11
        return BodyCenteredOrthorhombic()
    elseif ibrav == 12
        return PrimitiveMonoclinic()
    elseif ibrav == -12
        return PrimitiveMonoclinic()
    elseif ibrav == 13
        return CCenteredMonoclinic()
    elseif ibrav == -13
        return BCenteredMonoclinic()  # New in QE 6.5
    elseif ibrav == 14
        return PrimitiveTriclinic()
    else
        error("Bravais lattice undefined for `ibrav = $ibrav`!")
    end
end

"""
    Lattice(::Bravais, p[, obverse::Bool])

Create a Bravais lattice from the exact lattice type and cell parameters `p` (not `celldm`!).

The first elements of `p` are `a`, `b`, `c`; the last 3 are `α`, `β`, `γ` (in radians).
"""
Lattice(::PrimitiveCubic, p, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
])
Lattice(::FaceCenteredCubic, p, args...) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
])
function Lattice(::BodyCenteredCubic, p, obverse::Bool = true)
    if obverse
        return Lattice(p[1] / 2 * [
            1 1 1
            -1 1 1
            -1 -1 1
        ])
    else
        return Lattice(p[1] / 2 * [
            -1 1 1
            1 -1 1
            1 1 -1
        ])
    end
end
Lattice(::PrimitiveHexagonal, p::_Celldm, args...) = Lattice(p[1] * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 p[3]
])
function Lattice(::RCenteredHexagonal, p, obverse::Bool = true)
    if obverse
        r = p[4]
        tx = sqrt((1 - r) / 2)
        ty = sqrt((1 - r) / 6)
        tz = sqrt((1 + 2r) / 3)
        return Lattice(p[1] * [
            tx -ty tz
            0 2ty tz
            -tx -ty tz
        ])
    else  # -5
        ap = p[1] / √3
        γ = acos(p[4])
        ty = sqrt((1 - γ) / 6)
        tz = sqrt((1 + 2γ) / 3)
        u = tz - 2 * √2 * ty
        v = tz + √2 * ty
        return Lattice(ap * [
            u v v
            v u v
            v v u
        ])
    end
end
Lattice(::PrimitiveTetragonal, p::_Celldm, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3]
])
function Lattice(::BodyCenteredTetragonal, celldm::_Celldm, args...)
    r = celldm[3]
    return Lattice(celldm[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ])
end
Lattice(::PrimitiveOrthorhombic, p::_Celldm, args...) = Lattice(p[1] * [
    1 0 0
    0 p[2] 0
    0 0 p[3]
])
function Lattice(::BCenteredOrthorhombic, p::_Celldm, obverse::Bool = true)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    if obverse
        return Lattice([
            a / 2 b / 2 0
            -a / 2 b / 2 0
            0 0 c
        ])
    else
        return Lattice([
            a / 2 -b / 2 0
            a / 2 b / 2 0
            0 0 c
        ])
    end
end
Lattice(::ACenteredOrthorhombic, p, args...) = Lattice([
    p[1] 0 0
    0 p[2] / 2 -p[3] / 2
    0 p[2] / 2 p[3] / 2
])  # New in QE 6.4
function Lattice(::FaceCenteredOrthorhombic, p::_Celldm, args...)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2)
end
function Lattice(::BodyCenteredOrthorhombic, p::_Celldm, args...)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2)
end
# function Lattice(::PrimitiveMonoclinic, p, obverse::Bool = true)
#     a, b, c, _, β, γ = p[1:6]
#     if obverse
#         return Lattice([
#             a 0 0
#             b * cos(γ) b * sin(γ) 0
#             0 0 c
#         ])
#     else
#         return Lattice([
#             a 0 0
#             0 b 0
#             c * cos(β) 0 c * sin(β)
#         ])
#     end
# end
# function Lattice(::CCenteredMonoclinic, p, args...)
#     a, b, c = p[1:3]
#     return Lattice([
#         a / 2 0 -c / 2
#         b * cos(p[6]) b * sin(p[6]) 0
#         a / 2 0 c / 2
#     ])
# end
# function Lattice(::BCenteredMonoclinic, p, args...)
#     a, b, c = p[1:3]
#     return Lattice([
#         a / 2 b / 2 0
#         -a / 2 b / 2 0
#         c * cos(p[5]) 0 c * sin(p[5])
#     ])
# end
# function Lattice(::PrimitiveTriclinic, p, args...)
#     a, b, c, α, β, γ = p  # Every `p` that is an iterable can be used
#     δ = c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2) / sin(γ)
#     return Lattice([
#         a 0 0
#         b * cos(γ) b * sin(γ) 0
#         c * cos(β) c * (cos(α) - cos(β) * cos(γ)) / sin(γ) δ
#     ])
# end

end
