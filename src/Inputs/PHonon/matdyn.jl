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
