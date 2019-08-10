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

export ControlNamelist

# The following default values are picked from `<QE source>/Modules/read_namelists.f90`
@with_kw struct ControlNamelist <: Namelist
    calculation::String = "cp"
    title::String = "MD Simulation"
    verbosity::String = "low"
    isave::Int = 100
    restart_mode::String = "restart"
    wf_collect::Bool = true
    nstep::Int = 1
    iprint::Int = 10
    tstress::Bool = false
    tprnfor::Bool = false
    dt::Float64 = 1.0
    outdir::String = "./"
    wfcdir::String = "./"
    prefix::String = "cp"
    lkpoint_dir::Bool = true
    max_seconds::Float64 = 10000000.0
    etot_conv_thr::Float64 = 0.0001
    forc_conv_thr::Float64 = 0.001
    disk_io::String = "medium"
    pseudo_dir::String = "$(ENV["HOME"])/espresso/pseudo/"
    tefield::Bool = false
    dipfield::Bool = false
    lelfield::Bool = false
    nberrycyc::Int = 1
    lorbm::Bool = false
    lberry::Bool = false
    gdir::Int = 0
    nppstr::Int = 0
    lfcpopt::Bool = false
    gate::Bool = false
end  # struct ControlNamelist

end
