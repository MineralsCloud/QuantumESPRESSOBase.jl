const DVSCF_STAR = @NamedTuple begin
    open::Bool
    dir::String
    ext::String
    basis::String
    pat::Bool
end
const DRHO_STAR = @NamedTuple begin
    open::Bool
    dir::String
    ext::String
    basis::String
    pat::Bool
end

# The following default values are picked from `<QE source>/test-suite/not_epw_comp/phq_readin.f90`
"""
    PhNamelist <: Namelist

Represent the `INPUTPH` namelist of `ph.x`.
"""
@auto_hash_equals struct PhNamelist <: Namelist
    amass::Vector{Union{Nothing,Float64}}
    outdir::String
    prefix::String
    niter_ph::Int
    tr2_ph::Float64
    alpha_mix::Vector{Union{Nothing,Float64}}
    nmix_ph::Int
    verbosity::String
    reduce_io::Bool
    max_seconds::Float64
    fildyn::String
    fildrho::String
    fildvscf::String
    epsil::Bool
    lrpa::Bool
    lnoloc::Bool
    trans::Bool
    lraman::Bool
    eth_rps::Float64
    eth_ns::Float64
    dek::Float64
    recover::Bool
    low_directory_check::Bool
    only_init::Bool
    qplot::Bool
    q2d::Bool
    q_in_band_form::Bool
    electron_phonon::String
    lshift_q::Bool
    zeu::Bool  # The default value in QE's source code is `true`
    zue::Bool
    elop::Bool
    fpol::Bool
    ldisp::Bool
    nogg::Bool
    asr::Bool
    ldiag::Bool
    lqdir::Bool
    search_sym::Bool
    nq1::Int
    nq2::Int
    nq3::Int
    nk1::Int
    nk2::Int
    nk3::Int
    k1::Int
    k2::Int
    k3::Int
    start_irr::Int  # The default value in QE's source code is 0
    last_irr::Int
    nat_todo::Int
    modenum::Int
    start_q::Int
    last_q::Int
    dvscf_star::DVSCF_STAR
    drho_star::DRHO_STAR
end # struct PhNamelist
function PhNamelist(;
    amass = [0.0],
    outdir = "./",
    prefix = "pwscf",
    niter_ph = 100,
    tr2_ph = 1e-12,
    alpha_mix = 0.7 * ones(niter_ph),
    nmix_ph = 4,
    verbosity = "default",
    reduce_io = false,
    max_seconds = 1e7,
    fildyn = "matdyn",
    fildrho = " ",
    fildvscf = " ",
    epsil = false,
    lrpa = false,
    lnoloc = false,
    trans = true,
    lraman = false,
    eth_rps = 1e-9,
    eth_ns = 1e-12,
    dek = 1e-3,
    recover = false,
    low_directory_check = false,
    only_init = false,
    qplot = false,
    q2d = false,
    q_in_band_form = false,
    electron_phonon = " ",
    lshift_q = false,
    zeu = epsil,  # The default value in QE's source code is `true`
    zue = false,
    elop = false,
    fpol = false,
    ldisp = false,
    nogg = false,
    asr = false,
    ldiag = false,
    lqdir = false,
    search_sym = true,
    nq1 = 0,
    nq2 = 0,
    nq3 = 0,
    nk1 = 0,
    nk2 = 0,
    nk3 = 0,
    k1 = 0,
    k2 = 0,
    k3 = 0,
    start_irr = 1,  # The default value in QE's source code is 0
    last_irr = -1000,
    nat_todo = 0,
    modenum = 0,
    start_q = 1,
    last_q = -1000,
    dvscf_star = (
        open = false,
        dir = "./",
        ext = "dvscf",
        basis = "cartesian",
        pat = false,
    ),
    drho_star = (open = false, dir = "./", ext = "drho", basis = "modes", pat = true),
)
    return PhNamelist(
        amass,
        outdir,
        prefix,
        niter_ph,
        tr2_ph,
        alpha_mix,
        nmix_ph,
        verbosity,
        reduce_io,
        max_seconds,
        fildyn,
        fildrho,
        fildvscf,
        epsil,
        lrpa,
        lnoloc,
        trans,
        lraman,
        eth_rps,
        eth_ns,
        dek,
        recover,
        low_directory_check,
        only_init,
        qplot,
        q2d,
        q_in_band_form,
        electron_phonon,
        lshift_q,
        zeu,
        zue,
        elop,
        fpol,
        ldisp,
        nogg,
        asr,
        ldiag,
        lqdir,
        search_sym,
        nq1,
        nq2,
        nq3,
        nk1,
        nk2,
        nk3,
        k1,
        k2,
        k3,
        start_irr,
        last_irr,
        nat_todo,
        modenum,
        start_q,
        last_q,
        dvscf_star,
        drho_star,
    )
end
PhNamelist(nml::PhNamelist; kwargs...) = setproperties(nml, kwargs...)
PhNamelist(nml::PhNamelist, t::NamedTuple) = setproperties(nml, t)
PhNamelist(nml::PhNamelist, dict::AbstractDict) = setproperties(nml, dict)
