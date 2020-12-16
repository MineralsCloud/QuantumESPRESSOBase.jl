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
