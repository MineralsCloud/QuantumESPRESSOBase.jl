using AbInitioSoftwareBase.Inputs: Namelist
using ConstructionBase: setproperties

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
    PhNamelist(amass, outdir, prefix, niter_ph, tr2_ph, alpha_mix, nmix_ph, verbosity, reduce_io, max_seconds, fildyn, fildrho, fildvscf, epsil, lrpa, lnoloc, trans, lraman, eth_rps, eth_ns, dek, recover, low_directory_check, only_init, qplot, q2d, q_in_band_form, electron_phonon, lshift_q, zeu, zue, elop, fpol, ldisp, nogg, asr, ldiag, lqdir, search_sym, nq1, nq2, nq3, nk1, nk2, nk3, k1, k2, k3, start_irr, last_irr, nat_todo, modenum, start_q, last_q, dvscf_star, drho_star)
    PhNamelist(; kwargs...)
    PhNamelist(::PhNamelist; kwargs...)

Represent the `INPUTPH` namelist of `ph.x`.
"""
@struct_hash_equal struct PhNamelist <: Namelist
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
    amass=[0.0],
    outdir="./",
    prefix="pwscf",
    niter_ph=100,
    tr2_ph=1e-12,
    alpha_mix=0.7 * ones(niter_ph),
    nmix_ph=4,
    verbosity="default",
    reduce_io=false,
    max_seconds=1e7,
    fildyn="matdyn",
    fildrho=" ",
    fildvscf=" ",
    epsil=false,
    lrpa=false,
    lnoloc=false,
    trans=true,
    lraman=false,
    eth_rps=1e-9,
    eth_ns=1e-12,
    dek=1e-3,
    recover=false,
    low_directory_check=false,
    only_init=false,
    qplot=false,
    q2d=false,
    q_in_band_form=false,
    electron_phonon=" ",
    lshift_q=false,
    zeu=epsil,  # The default value in QE's source code is `true`
    zue=false,
    elop=false,
    fpol=false,
    ldisp=false,
    nogg=false,
    asr=false,
    ldiag=false,
    lqdir=false,
    search_sym=true,
    nq1=0,
    nq2=0,
    nq3=0,
    nk1=0,
    nk2=0,
    nk3=0,
    k1=0,
    k2=0,
    k3=0,
    start_irr=1,  # The default value in QE's source code is 0
    last_irr=-1000,
    nat_todo=0,
    modenum=0,
    start_q=1,
    last_q=-1000,
    dvscf_star=(open=false, dir="./", ext="dvscf", basis="cartesian", pat=false),
    drho_star=(open=false, dir="./", ext="drho", basis="modes", pat=true),
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

function (x::VerbositySetter)(control::PhNamelist)
    @set! control.verbosity = x.v
    return control
end

# The following default values are picked from `<QE source>/PHonon/PH/q2r.f90`
"""
    Q2rNamelist(fildyn, flfrc, loto_2d, zasr)
    Q2rNamelist(; kwargs...)
    Q2rNamelist(::Q2rNamelist; kwargs...)

Represent the `INPUT` namelist of `q2r.x`.
"""
@struct_hash_equal struct Q2rNamelist <: Namelist
    fildyn::String
    flfrc::String
    loto_2d::Bool
    zasr::String
end # struct Q2rNamelist
function Q2rNamelist(; fildyn=" ", flfrc=" ", loto_2d=false, zasr="no")
    return Q2rNamelist(fildyn, flfrc, loto_2d, zasr)
end
Q2rNamelist(nml::Q2rNamelist; kwargs...) = setproperties(nml, kwargs...)

# The following default values are picked from `<QE source>/PHonon/PH/matdyn.f90`
"""
    MatdynNamelist(dos, deltaE, ndos, nk1, nk2, nk3, asr, readtau, flfrc, fldos, flfrq, flvec, fleig, fldyn, fltau, amass, at, ntyp, l1, l2, l3, la2F, q_in_band_form, eigen_similarity, q_in_cryst_coord, na_ifc, fd, nosym, loto_2d)
    MatdynNamelist(; kwargs...)
    MatdynNamelist(::MatdynNamelist; kwargs...)

Represent the `INPUT` namelist of `matdyn.x`.
"""
@struct_hash_equal struct MatdynNamelist <: Namelist
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
    dos=false,
    deltaE=1.0,
    ndos=1,
    nk1=0,
    nk2=0,
    nk3=0,
    asr="no",
    readtau=false,
    flfrc=" ",
    fldos="matdyn.dos",
    flfrq="matdyn.freq",
    flvec="matdyn.modes",
    fleig="matdyn.eig",
    fldyn=" ",
    fltau=" ",
    amass=zeros(1),
    at=zeros(3, 3),  # FIXME: not very sure
    ntyp=0,
    l1=1,
    l2=1,
    l3=1,
    la2F=false,
    q_in_band_form=false,
    eigen_similarity=false,
    q_in_cryst_coord=false,
    na_ifc=false,
    fd=false,
    nosym=false,
    loto_2d=false,
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

"""
    DynmatNamelist(asr, axis, fildyn, filout, filmol, filxsf, fileig, amass, q, lperm, lplasma)
    DynmatNamelist(; kwargs...)
    DynmatNamelist(::DynmatNamelist; kwargs...)

Represent the `INPUT` namelist of `dynmat.x`.
"""
@struct_hash_equal struct DynmatNamelist <: Namelist
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
    asr="no",
    axis=3,
    fildyn="matdyn",
    filout="dynmat.out",
    filmol="dynmat.mold",
    filxsf="dynmat.axsf",
    fileig=" ",
    amass=zeros(1),
    q=zeros(3),
    lperm=false,
    lplasma=false,
)
    return DynmatNamelist(
        asr, axis, fildyn, filout, filmol, filxsf, fileig, amass, q, lperm, lplasma
    )
end
DynmatNamelist(nml::DynmatNamelist; kwargs...) = setproperties(nml, kwargs...)

groupname(::Type{PhNamelist}) = "INPUTPH"
groupname(::Type{Q2rNamelist}) = "INPUT"
groupname(::Type{MatdynNamelist}) = "INPUT"
groupname(::Type{DynmatNamelist}) = "INPUT"
