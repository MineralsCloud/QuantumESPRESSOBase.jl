"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using AbInitioSoftwareBase.Inputs: Namelist
using AutoHashEquals: @auto_hash_equals
using Compat: @NamedTuple
using ConstructionBase: setproperties
using Setfield: @set!

using ..Inputs: Card, QuantumESPRESSOInput
using ..Inputs.PWscf: SpecialPoint, PWInput

import ..Inputs

export SpecialPoint,
    QPointsCard,
    PhInput,
    Q2rInput,
    MatdynInput,
    DynmatInput,
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist
export relayinfo

include("ph.jl")
include("q2r.jl")
include("matdyn.jl")

"""
    DynmatNamelist <: Namelist

Represent the `INPUT` namelist of `dynmat.x`.
"""
@auto_hash_equals struct DynmatNamelist <: Namelist
    asr::String
    axis::Int
    fildyn::String
    filout::String
    filmol::String
    filxsf::String
    fileig::String
    amass::Vector{Union{Nothing,Float64}}
    q::Vector{Union{Nothing,Float64}}
    lperm::Bool
    lplasma::Bool
end # struct DynmatNamelist
function DynmatNamelist(;
    asr = "no",
    axis = 3,
    fildyn = "matdyn",
    filout = "dynmat.out",
    filmol = "dynmat.mold",
    filxsf = "dynmat.axsf",
    fileig = " ",
    amass = zeros(1),
    q = zeros(3),
    lperm = false,
    lplasma = false,
)
    return DynmatNamelist(
        asr,
        axis,
        fildyn,
        filout,
        filmol,
        filxsf,
        fileig,
        amass,
        q,
        lperm,
        lplasma,
    )
end
DynmatNamelist(nml::DynmatNamelist; kwargs...) = setproperties(nml, kwargs...)
DynmatNamelist(nml::DynmatNamelist, t::NamedTuple) = setproperties(nml, t)
DynmatNamelist(nml::DynmatNamelist, dict::AbstractDict) = setproperties(nml, dict)

Inputs.groupname(::Type{PhNamelist}) = "INPUTPH"
Inputs.groupname(::Type{Q2rNamelist}) = "INPUT"
Inputs.groupname(::Type{MatdynNamelist}) = "INPUT"
Inputs.groupname(::Type{DynmatNamelist}) = "INPUT"

struct QPointsCard <: Card
    data::Vector{SpecialPoint}
end

struct PhInput <: QuantumESPRESSOInput
    title_line::String
    inputph::PhNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct PhInput
PhInput(inputph::PhNamelist, qpts::QPointsCard) = PhInput(inputph.prefix, inputph, qpts)
PhInput(inputph::PhNamelist) = PhInput(inputph.prefix, inputph, nothing)
PhInput() = PhInput(PhNamelist().prefix, PhNamelist(), nothing)

struct Q2rInput <: QuantumESPRESSOInput
    input::Q2rNamelist
end # struct Q2rInput
Q2rInput() = Q2rInput(Q2rNamelist())

struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct MatdynInput
MatdynInput(input) = MatdynInput(input, nothing)
MatdynInput() = MatdynInput(MatdynNamelist(), nothing)

struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist
end # struct DynmatInput
DynmatInput() = DynmatInput(DynmatNamelist())

"""
    relayinfo(from::PWInput, to::PhInput)

Relay shared information from a `PWInput` to a `PhInput`.

A `PWInput` before a `PhInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function relayinfo(pw::PWInput, ph::PhInput)
    @set! ph.inputph.outdir = pw.control.outdir
    @set! ph.inputph.prefix = pw.control.prefix
    return ph
end # function relayinfo
"""
    relayinfo(from::PhInput, to::Q2rInput)

Relay shared information from a `PhInput` to a `Q2rInput`.

A `PhInput` before a `Q2rInput` has the information of `fildyn`. It must keep the same in a q2r calculation.
"""
function relayinfo(ph::PhInput, q2r::Q2rInput)
    @set! q2r.input.fildyn = ph.inputph.fildyn
    return q2r
end # function relayinfo
"""
    relayinfo(from::Q2rInput, to::MatdynInput)

Relay shared information from a `Q2rInput` to a `MatdynInput`.

A `Q2rInput` before a `MatdynInput` has the information of `fildyn`, `flfrc` and `loto_2d`. They must keep the same
in a matdyn calculation.
"""
function relayinfo(q2r::Q2rInput, matdyn::MatdynInput)
    @set! matdyn.input.flfrc = q2r.input.flfrc
    @set! matdyn.input.loto_2d = q2r.input.loto_2d
    return matdyn
end # function relayinfo
function relayinfo(ph::PhInput, matdyn::MatdynInput)
    @set! matdyn.input.amass = ph.inputph.amass
    @set! matdyn.input.q_in_band_form = ph.inputph.q_in_band_form
    return matdyn
end # function relayinfo
"""
    relayinfo(from::PhInput, to::DynmatInput)

Relay shared information from a `PhInput` to a `DynmatInput`.

A `PhInput` before a `DynmatInput` has the information of `asr`, `fildyn` and `amass`. They must keep the same
in a dynmat calculation.
"""
function relayinfo(ph::PhInput, dynmat::DynmatInput)
    # @set! dynmat.input.asr = ph.inputph.asr  # TODO
    @set! dynmat.input.fildyn = ph.inputph.fildyn
    @set! dynmat.input.amass = ph.inputph.amass
    return dynmat
end # function relayinfo

end
