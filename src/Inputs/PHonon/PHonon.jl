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

# The following default values are picked from `<QE source>/PHonon/PH/q2r.f90`
"""
    Q2rNamelist <: Namelist

Represent the `INPUT` namelist of `q2r.x`.
"""
@auto_hash_equals struct Q2rNamelist <: Namelist
    fildyn::String
    flfrc::String
    loto_2d::Bool
    zasr::String
end # struct Q2rNamelist
function Q2rNamelist(; fildyn = " ", flfrc = " ", loto_2d = false, zasr = "no")
    return Q2rNamelist(fildyn, flfrc, loto_2d, zasr)
end
Q2rNamelist(nml::Q2rNamelist; kwargs...) = setproperties(nml, kwargs...)
Q2rNamelist(nml::Q2rNamelist, t::NamedTuple) = setproperties(nml, t)
Q2rNamelist(nml::Q2rNamelist, dict::AbstractDict) = setproperties(nml, dict)

# The following default values are picked from `<QE source>/PHonon/PH/matdyn.f90`
"""
    MatdynNamelist <: Namelist

Represent the `INPUT` namelist of `matdyn.x`.
"""
@auto_hash_equals struct MatdynNamelist <: Namelist
    dos::Bool
    deltaE::Float64
    ndos::Int
    nk1::Int
    nk2::Int
    nk3::Int
    asr::String
    readtau::Bool
    flfrc::String
    fldos::String
    flfrq::String
    flvec::String
    fleig::String
    fldyn::String
    fltau::String
    amass::Vector{Union{Nothing,Float64}}
    at::Matrix{Union{Nothing,Float64}}  # FIXME: not very sure
    ntyp::Int
    l1::Int
    l2::Int
    l3::Int
    la2F::Bool
    q_in_band_form::Bool
    eigen_similarity::Bool
    q_in_cryst_coord::Bool
    na_ifc::Bool
    fd::Bool
    nosym::Bool
    loto_2d::Bool
end # struct MatdynNamelist
function MatdynNamelist(;
    dos = false,
    deltaE = 1.0,
    ndos = 1,
    nk1 = 0,
    nk2 = 0,
    nk3 = 0,
    asr = "no",
    readtau = false,
    flfrc = " ",
    fldos = "matdyn.dos",
    flfrq = "matdyn.freq",
    flvec = "matdyn.modes",
    fleig = "matdyn.eig",
    fldyn = " ",
    fltau = " ",
    amass = zeros(1),
    at = zeros(3, 3),  # FIXME: not very sure
    ntyp = 0,
    l1 = 1,
    l2 = 1,
    l3 = 1,
    la2F = false,
    q_in_band_form = false,
    eigen_similarity = false,
    q_in_cryst_coord = false,
    na_ifc = false,
    fd = false,
    nosym = false,
    loto_2d = false,
)
    return MatdynNamelist(
        dos,
        deltaE,
        ndos,
        nk1,
        nk2,
        nk3,
        asr,
        readtau,
        flfrc,
        fldos,
        flfrq,
        flvec,
        fleig,
        fldyn,
        fltau,
        amass,
        at,
        ntyp,
        l1,
        l2,
        l3,
        la2F,
        q_in_band_form,
        eigen_similarity,
        q_in_cryst_coord,
        na_ifc,
        fd,
        nosym,
        loto_2d,
    )
end
MatdynNamelist(nml::MatdynNamelist; kwargs...) = setproperties(nml, kwargs...)
MatdynNamelist(nml::MatdynNamelist, t::NamedTuple) = setproperties(nml, t)
MatdynNamelist(nml::MatdynNamelist, dict::AbstractDict) = setproperties(nml, dict)

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
