"""
# module PHonon



# Examples

```jldoctest
julia>
```
"""
module PHonon

using Parameters: @with_kw

using ..Inputs: Namelist, Card, Input
using ..Inputs.PWscf: SpecialKPoint

import ..Inputs

export SpecialKPoint, QPointsCard, PhInput, Q2rInput, MatdynInput, DynmatInput
export PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist

# The following default values are picked from `<QE source>/test-suite/not_epw_comp/phq_readin.f90`
"""
    PhNamelist <: Namelist

Represent the `INPUTPH` namelist of `ph.x`.
"""
@with_kw struct PhNamelist <: Namelist
    amass::Vector{Union{Nothing,Float64}} = [0.0]
    outdir::String = "./"
    prefix::String = "pwscf"
    niter_ph::Int = 100
    tr2_ph::Float64 = 1e-12
    alpha_mix::Vector{Union{Nothing,Float64}} = 0.7 * ones(niter_ph)
    nmix_ph::Int = 4
    verbosity::String = "default"
    reduce_io::Bool = false
    max_seconds::Float64 = 1e7
    fildyn::String = "matdyn"
    fildrho::String = " "
    fildvscf::String = " "
    epsil::Bool = false
    lrpa::Bool = false
    lnoloc::Bool = false
    trans::Bool = true
    lraman::Bool = false
    eth_rps::Float64 = 1e-9
    eth_ns::Float64 = 1e-12
    dek::Float64 = 1e-3
    recover::Bool = false
    low_directory_check::Bool = false
    only_init::Bool = false
    qplot::Bool = false
    q2d::Bool = false
    q_in_band_form::Bool = false
    electron_phonon::String = " "
    lshift_q::Bool = false
    zeu::Bool = epsil  # The default value in QE's source code is `true`
    zue::Bool = false
    elop::Bool = false
    fpol::Bool = false
    ldisp::Bool = false
    nogg::Bool = false
    asr::Bool = false
    ldiag::Bool = false
    lqdir::Bool = false
    search_sym::Bool = true
    nq1::Int = 0
    nq2::Int = 0
    nq3::Int = 0
    nk1::Int = 0
    nk2::Int = 0
    nk3::Int = 0
    k1::Int = 0
    k2::Int = 0
    k3::Int = 0
    start_irr::Int = 1  # The default value in QE's source code is 0
    last_irr::Int = -1000
    nat_todo::Int = 0
    modenum::Int = 0
    start_q::Int = 1
    last_q::Int = -1000
    # dvscf_star::String = 1
    # drho_star::String = 1
end # struct PhNamelist

# The following default values are picked from `<QE source>/PHonon/PH/q2r.f90`
"""
    Q2rNamelist <: Namelist

Represent the `INPUT` namelist of `q2r.x`.
"""
@with_kw struct Q2rNamelist <: Namelist
    fildyn::String = " "
    flfrc::String = " "
    loto_2d::Bool = false
    zasr::String = "no"
end # struct Q2rNamelist

# The following default values are picked from `<QE source>/PHonon/PH/matdyn.f90`
"""
    MatdynNamelist <: Namelist

Represent the `INPUT` namelist of `matdyn.x`.
"""
@with_kw struct MatdynNamelist <: Namelist
    dos::Bool = false
    deltaE::Float64 = 1.0
    ndos::Int = 1
    nk1::Int = 0
    nk2::Int = 0
    nk3::Int = 0
    asr::String = "no"
    readtau::Bool = false
    flfrc::String = " "
    fldos::String = "matdyn.dos"
    flfrq::String = "matdyn.freq"
    flvec::String = "matdyn.modes"
    fleig::String = "matdyn.eig"
    fldyn::String = " "
    fltau::String = " "
    amass::Vector{Union{Nothing,Float64}} = zeros(1)
    at::Matrix{Union{Nothing,Float64}} = zeros(3, 3)  # FIXME: not very sure
    ntyp::Int = 0
    l1::Int = 1
    l2::Int = 1
    l3::Int = 1
    la2F::Bool = false
    q_in_band_form::Bool = false
    eigen_similarity::Bool = false
    q_in_cryst_coord::Bool = false
    na_ifc::Bool = false
    fd::Bool = false
    nosym::Bool = false
    loto_2d::Bool = false
end # struct MatdynNamelist

"""
    DynmatNamelist <: Namelist

Represent the `INPUT` namelist of `dynmat.x`.
"""
@with_kw struct DynmatNamelist <: Namelist
    asr::String = "no"
    axis::Int = 3
    fildyn::String = "matdyn"
    filout::String = "dynmat.out"
    filmol::String = "dynmat.mold"
    filxsf::String = "dynmat.axsf"
    fileig::String = " "
    amass::Vector{Union{Nothing,Float64}} = zeros(1)
    q::Vector{Union{Nothing,Float64}} = zeros(3)
    lperm::Bool = false
    lplasma::Bool = false
end # struct DynmatNamelist

Inputs.titleof(::Type{PhNamelist}) = "INPUTPH"
Inputs.titleof(::Type{Q2rNamelist}) = "INPUT"
Inputs.titleof(::Type{MatdynNamelist}) = "INPUT"
Inputs.titleof(::Type{DynmatNamelist}) = "INPUT"

struct QPointsCard <: Card
    data::Vector{SpecialKPoint}
end

struct PhInput <: Input
    inputph::PhNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct PhInput
PhInput(inputph = PhNamelist(), q_points = nothing) = PhInput(inputph, q_points)

struct Q2rInput <: Input
    input::Q2rNamelist
end # struct Q2rInput
Q2rInput() = Q2rInput(Q2rNamelist())

struct MatdynInput <: Input
    input::MatdynNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct MatdynInput
MatdynInput(input = MatdynNamelist(), q_points = nothing) = MatdynInput(input, q_points)

struct DynmatInput <: Input
    input::DynmatNamelist
end # struct DynmatInput
DynmatInput() = DynmatInput(DynmatNamelist())

end
