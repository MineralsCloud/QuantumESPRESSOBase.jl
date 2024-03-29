using ..QuantumESPRESSOBase: QuantumESPRESSOInput
using ..PWscf: PWInput

@struct_hash_equal struct PhInput <: QuantumESPRESSOInput
    title_line::String
    inputph::PhNamelist
    q_points::Union{Nothing,QPointsCard}
end
PhInput(inputph::PhNamelist, qpts::QPointsCard) = PhInput(inputph.prefix, inputph, qpts)
PhInput(inputph::PhNamelist) = PhInput(inputph.prefix, inputph, nothing)
PhInput() = PhInput(PhNamelist().prefix, PhNamelist(), nothing)

struct Q2rInput <: QuantumESPRESSOInput
    input::Q2rNamelist
end
Q2rInput() = Q2rInput(Q2rNamelist())

@struct_hash_equal struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist
    q_points::Union{Nothing,QPointsCard}
end
MatdynInput(input) = MatdynInput(input, nothing)
MatdynInput() = MatdynInput(MatdynNamelist(), nothing)

@struct_hash_equal struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist
end
DynmatInput() = DynmatInput(DynmatNamelist())

function (x::VerbositySetter)(template::PhInput)
    @reset template.inputph.verbosity = x.v
    return template
end

"""
    relayinfo(from::PWInput, to::PhInput)

Relay shared information from a `PWInput` to a `PhInput`.

A `PWInput` before a `PhInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function relayinfo(pw::PWInput, ph::PhInput)
    @reset ph.inputph.outdir = pw.control.outdir
    @reset ph.inputph.prefix = pw.control.prefix
    return ph
end
"""
    relayinfo(from::PhInput, to::Q2rInput)

Relay shared information from a `PhInput` to a `Q2rInput`.

A `PhInput` before a `Q2rInput` has the information of `fildyn`. It must keep the same in a q2r calculation.
"""
function relayinfo(ph::PhInput, q2r::Q2rInput)
    @reset q2r.input.fildyn = ph.inputph.fildyn
    return q2r
end
"""
    relayinfo(from::Q2rInput, to::MatdynInput)

Relay shared information from a `Q2rInput` to a `MatdynInput`.

A `Q2rInput` before a `MatdynInput` has the information of `fildyn`, `flfrc` and `loto_2d`. They must keep the same
in a matdyn calculation.
"""
function relayinfo(q2r::Q2rInput, matdyn::MatdynInput)
    @reset matdyn.input.flfrc = q2r.input.flfrc
    @reset matdyn.input.loto_2d = q2r.input.loto_2d
    return matdyn
end
function relayinfo(ph::PhInput, matdyn::MatdynInput)
    @reset matdyn.input.amass = ph.inputph.amass
    @reset matdyn.input.q_in_band_form = ph.inputph.q_in_band_form
    return matdyn
end
"""
    relayinfo(from::PhInput, to::DynmatInput)

Relay shared information from a `PhInput` to a `DynmatInput`.

A `PhInput` before a `DynmatInput` has the information of `asr`, `fildyn` and `amass`. They must keep the same
in a dynmat calculation.
"""
function relayinfo(ph::PhInput, dynmat::DynmatInput)
    # @reset dynmat.input.asr = ph.inputph.asr  # TODO
    @reset dynmat.input.fildyn = ph.inputph.fildyn
    @reset dynmat.input.amass = ph.inputph.amass
    return dynmat
end
