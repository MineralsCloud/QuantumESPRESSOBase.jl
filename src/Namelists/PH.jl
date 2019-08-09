"""
# module PH



# Examples

```jldoctest
julia>
```
"""
module PH

using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists: Namelist

export INPUTPHNamelist

@with_kw struct INPUTPHNamelist <: Namelist
    amass::Vector{Union{Missing, Float64}} = [0.0]
    outdir::String = "./"
    prefix::String = "pwscf"
    niter_ph::Int = 100
    tr2_ph::Float64 = 1e-12
    alpha_mix::Vector{Union{Missing, Float64}} = 0.7 * ones(niter_ph)
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
end  # struct INPUTPHNamelist

end
