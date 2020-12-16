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
