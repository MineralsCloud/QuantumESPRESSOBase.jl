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
using CrystallographyBase: ReciprocalPoint
using Setfield: @set!

using ..Inputs: Card, QuantumESPRESSOInput, VerbositySetter
using ..Inputs.PWscf: PWInput

import ..Inputs: groupname

export QPointsCard,
    PhInput,
    Q2rInput,
    MatdynInput,
    DynmatInput,
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist,
    VerbositySetter
export relayinfo

include("ph.jl")
include("q2r.jl")
include("matdyn.jl")
include("dynmat.jl")

groupname(::Type{PhNamelist}) = "INPUTPH"
groupname(::Type{Q2rNamelist}) = "INPUT"
groupname(::Type{MatdynNamelist}) = "INPUT"
groupname(::Type{DynmatNamelist}) = "INPUT"

struct QPointsCard <: Card
    data::Vector{ReciprocalPoint}
end

include("inputs.jl")

function (x::VerbositySetter)(control::PhNamelist)
    @set! control.verbosity = x.v
    return control
end
function (x::VerbositySetter)(template::PhInput)
    @set! template.inputph.verbosity = x.v
    return template
end

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
