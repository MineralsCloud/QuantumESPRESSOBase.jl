"""
# module CP



# Examples

```jldoctest
julia>
```
"""
module CP

using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists: Namelist

export ControlNamelist, SystemNamelist

# The following default values are picked from `<QE source>/Modules/read_namelists.f90`
@with_kw struct ControlNamelist <: Namelist
    calculation::String = "cp"
    title::String = "MD Simulation"
    verbosity::String = "low"
    isave::Int = 100
    restart_mode::String = "restart"
    nstep::Int = 50
    iprint::Int = 10
    tstress::Bool = false
    tprnfor::Bool = false
    dt::Float64 = 1.0
    outdir::String = "./"
    saverho::Bool = true
    prefix::String = "cp"
    ndr::Int = 50
    ndw::Int = 50
    tabps::Bool = false
    max_seconds::Float64 = 1e7
    etot_conv_thr::Float64 = 1e-4
    forc_conv_thr::Float64 = 1e-3
    ekin_conv_thr::Float64 = 1e-6
    disk_io::String = "default"
    memory::String = "default"
    pseudo_dir::String = "$(ENV["HOME"])/espresso/pseudo/"
    tefield::Bool = false
end  # struct ControlNamelist

@with_kw struct SystemNamelist <: Namelist
    ibrav::Int = -1
    celldm::Vector{Union{Missing, Float64}} = zeros(6); @assert length(celldm) ≤ 6
    A::Float64 = 0.0
    B::Float64 = 0.0
    C::Float64 = 0.0
    cosAB::Float64 = 0.0
    cosAC::Float64 = 0.0
    cosBC::Float64 = 0.0
    nat::Int = 0
    ntyp::Int = 0
    nbnd::Int = 0
    tot_charge::Float64 = 0.0
    tot_magnetization::Float64 = -1.0
    ecutwfc::Float64 = 0.0
    ecutrho::Float64 = 0.0
    nr1::Int = 0
    nr2::Int = 0
    nr3::Int = 0
    nr1s::Int = 0
    nr2s::Int = 0
    nr3s::Int = 0
    nr1b::Int = 0
    nr2b::Int = 0
    nr3b::Int = 0
    occupations::String = "fixed"
    degauss::Float64 = 0.0
    smearing::String = "gaussian"
    nspin::Int = 1
    ecfixed::Float64 = 0.0
    qcutz::Float64 = 0.0
    q2sigma::Float64 = 0.1  # The default value in QE's source code is 0.01
    input_dft::String = "none"
    exx_fraction::Float64 = 0.25
    lda_plus_u::Bool = false
    Hubbard_U::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    vdw_corr::String = "none"
    london_s6::Float64 = 0.75
    london_rcut::Float64 = 200.0
    ts_vdw::Bool = false
    ts_vdw_econv_thr::Float64 = 1e-6
    ts_vdw_isolated::Bool = false
    assume_isolated::Bool = "none"
end  # struct SystemNamelist

end
