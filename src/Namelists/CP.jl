"""
# module CP



# Examples

```jldoctest
julia>
```
"""
module CP

using LinearAlgebra: det

using ConstructionBase: setproperties
using Parameters: @with_kw

using QuantumESPRESSOBase.Setters: VerbositySetter
using QuantumESPRESSOBase.Namelists: Namelist

import QuantumESPRESSOBase
import QuantumESPRESSOBase.Setters

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    PressAiNamelist,
    WannierNamelist

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
    pseudo_dir::String = raw"$HOME/espresso/pseudo/"
    tefield::Bool = false
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1282-L1369.
    @assert(
        calculation ∈ (
            "cp",
            "scf",
            "nscf",
            "relax",
            "vc-relax",
            "vc-cp",
            "cp-wf",
            "vc-cp-wf",
        )
    )
    @assert verbosity ∈ ("high", "low", "debug", "medium", "default", "minimal")
    @assert isave >= 1
    @assert nstep >= 0
    @assert iprint >= 1
    @assert dt >= 0
    @assert ndr >= 50
    @assert ndw <= 0 || ndw >= 50
    @assert max_seconds >= 0
    @assert etot_conv_thr >= 0
    @assert forc_conv_thr >= 0
    @assert ekin_conv_thr >= 0
    @assert memory ∈ ("small", "default", "large")
end # struct ControlNamelist

@with_kw struct SystemNamelist <: Namelist
    ibrav::Int = -1
    celldm::Vector{Union{Nothing,Float64}} = Vector{Float64}(undef, 6)
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
    Hubbard_U::Vector{Union{Nothing,Float64}} = []
    vdw_corr::String = "none"
    london_s6::Float64 = 0.75
    london_rcut::Float64 = 200.0
    ts_vdw::Bool = false
    ts_vdw_econv_thr::Float64 = 1e-6
    ts_vdw_isolated::Bool = false
    assume_isolated::String = "none"
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1378-L1499.
    @assert(
        if ibrav == 0
            celldm[1] != 0 || A != 0  # Cannot use `iszero` to compare now!
        else
            true
        end,
        "Invalid lattice parameters (`celldm` $celldm or `A` $A)!"
    )
    @assert(
        if ibrav == 14
            length(celldm) == 6
        elseif ibrav ∈ (5, -5, 12, 13)
            4 <= length(celldm) <= 6
        elseif ibrav ∈ (4, 6, 7, 8, 9, -9, 10, 11)
            3 <= length(celldm) <= 6
        else
            1 <= length(celldm) <= 6
        end,
        "`celldm` must have length between 1 to 6! See `ibrav`'s doc!"
    )
    @assert nat >= 0
    @assert(0 <= ntyp <= 10, "`ntyp` $ntyp is either less than zero or too large!")
    @assert nspin ∈ (1, 2, 4)
    @assert ecutwfc >= 0
    @assert ecutrho >= 0
    @assert(iszero(degauss), "`degauss` is not used in CP!")
    @assert ecfixed >= 0
    @assert qcutz >= 0
    @assert q2sigma >= 0
end # struct SystemNamelist

@with_kw struct ElectronsNamelist <: Namelist
    electron_maxstep::Int = 100
    electron_dynamics::String = "none"
    conv_thr::Float64 = 1e-6
    niter_cg_restart::Int = 20
    efield::Float64 = 0.0
    epol::Int = 3
    emass::Float64 = 400.0
    emass_cutoff::Float64 = 2.5
    orthogonalization::String = "ortho"
    ortho_eps::Float64 = 1e-8
    ortho_max::Int = 20
    ortho_para::Int = 0
    electron_damping::Float64 = 0.1
    electron_velocities::String = "default"
    electron_temperature::String = "not_controlled"
    ekincw::Float64 = 0.001
    fnosee::Float64 = 1.0
    startingwfc::String = "random"
    tcg::Bool = false
    maxiter::Int = 100
    passop::Float64 = 0.3
    n_inner::Int = 2
    ninter_cold_restart::Int = 1
    lambda_cold::Float64 = 0.03
    grease::Float64 = 1.0
    ampre::Float64 = 0.0
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1508-L1543.
    @assert electron_dynamics ∈ ("none", "sd", "damp", "verlet", "cg")  # Different from code
    @assert emass > 0
    @assert emass_cutoff > 0
    @assert orthogonalization ∈ ("ortho", "Gram-Schmidt")
    @assert ortho_eps > 0
    @assert ortho_max >= 1
    @assert electron_velocities ∈ ("zero", "default", "change_step")  # New in 6.4
    @assert electron_temperature ∈ ("nose", "rescaling", "not_controlled")
    @assert ekincw > 0
    @assert fnosee > 0
    @assert startingwfc ∈ ("atomic", "random")
end # struct ElectronsNamelist

@with_kw struct IonsNamelist <: Namelist
    ion_dynamics::String = "none"
    ion_positions::String = "default"
    ion_velocities::String = "default"
    ion_damping::Float64 = 0.2
    ion_radius::Vector{Union{Nothing,Float64}} = [0.5]  # The default value in QE's source code is just one 0.5
    iesr::Int = 1
    ion_nstepe::Int = 1
    remove_rigid_rot::Bool = false
    ion_temperature::String = "not_controlled"
    tempw::Float64 = 300.0
    fnosep::Float64 = 1.0
    tolp::Float64 = 100.0
    nhpcl::Int = 1
    nhptyp::Int = 0
    nhgrp::Vector{Union{Nothing,Int}} = [0]
    fnhscl::Vector{Union{Nothing,Float64}} = zeros(1)
    ndega::Int = 0
    tranp::Vector{Union{Nothing,Bool}} = falses(1)
    amprp::Vector{Union{Nothing,Float64}} = zeros(1)
    greasp::Float64 = 1.0
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1552-L1585.
    @assert ion_dynamics ∈ ("none", "sd", "cg", "damp", "verlet")
    @assert ion_positions ∈ ("default", "from_input")
    @assert ion_velocities ∈ ("default", "change_step", "random", "from_input", "zero")
    @assert ion_nstepe > 0
    @assert ion_temperature ∈ ("nose", "rescaling", "not_controlled")
    @assert tempw > 0
    @assert fnosep > 0
    @assert 0 <= nhpcl <= 4

end # struct IonsNamelist

@with_kw struct CellNamelist <: Namelist
    cell_parameters::String = "default"
    cell_dynamics::String = "none"
    cell_velocities::String = "default"
    cell_damping::Float64 = 0.1
    press::Float64 = 0.0
    wmass::Float64 = 0.0
    cell_factor::Float64 = 1.2
    cell_temperature::String = "not_controlled"
    temph::Float64 = 0.0
    fnoseh::Float64 = 1.0
    greash::Float64 = 1.0
    cell_dofree::String = "all"
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1596-L1625.
    @assert cell_parameters ∈ ("default", "from_input")
    @assert cell_dynamics ∈ ("none", "sd", "damp-pr", "pr")
    @assert cell_velocities ∈ ("zero", "default")
    @assert wmass >= 0
    @assert cell_temperature ∈ ("nose", "rescaling", "not_controlled")
    @assert(
        cell_dofree ∈ (
            "all",
            "x",
            "y",
            "z",
            "xy",
            "xz",
            "yz",
            "xyz",
            "shape",
            "2Dxy",
            "2Dshape",
            "volume",
        )
    )
end # struct CellNamelist

@with_kw struct PressAiNamelist <: Namelist
    abivol::Bool = false
    abisur::Bool = false
    P_ext::Float64 = 0.0
    pvar::Bool = false
    P_in::Float64 = 0.0
    P_fin::Float64 = 0.0
    Surf_t::Float64 = 0.0
    rho_thr::Float64 = 0.0
    dthr::Float64 = 0.0
end # struct PressAiNamelist

@with_kw struct WannierNamelist <: Namelist
    wf_efield::Bool = false
    wf_switch::Bool = false
    sw_len::Int = 1
    efx0::Float64 = 0.0
    efy0::Float64 = 0.0
    efz0::Float64 = 0.0
    efx1::Float64 = 0.0
    efy1::Float64 = 0.0
    efz1::Float64 = 0.0
    wfsd::Int = 1
    wfdt::Float64 = 5.0
    maxwfdt::Float64 = 0.3
    nit::Int = 10
    nsd::Int = 10
    wf_q::Float64 = 1500.0
    wf_friction::Float64 = 0.3
    nsteps::Int = 20
    tolw::Float64 = 1e-8
    adapt::Bool = true
    calwf::Int = 3
    nwf::Int = 0
    wffort::Int = 40
    writev::Bool = false
    exx_neigh::Int = 60
    exx_dis_cutoff::Float64 = 8.0
    exx_poisson_eps::Float64 = 1e-6
    exx_ps_rcut_self::Float64 = 6.0
    exx_ps_rcut_pair::Float64 = 5.0
    exx_me_rcut_self::Float64 = 10.0
    exx_me_rcut_pair::Float64 = 7.0
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1634-L1650.
    @assert 1 <= calwf <= 5
    @assert 1 <= wfsd <= 3
end # struct WannierNamelist

QuantumESPRESSOBase.asfieldname(::Type{<:ControlNamelist}) = :control
QuantumESPRESSOBase.asfieldname(::Type{<:SystemNamelist}) = :system
QuantumESPRESSOBase.asfieldname(::Type{<:ElectronsNamelist}) = :electrons
QuantumESPRESSOBase.asfieldname(::Type{<:IonsNamelist}) = :ions
QuantumESPRESSOBase.asfieldname(::Type{<:CellNamelist}) = :cell

QuantumESPRESSOBase.titleof(::Type{<:ControlNamelist}) = "CONTROL"
QuantumESPRESSOBase.titleof(::Type{<:SystemNamelist}) = "SYSTEM"
QuantumESPRESSOBase.titleof(::Type{<:ElectronsNamelist}) = "ELECTRONS"
QuantumESPRESSOBase.titleof(::Type{<:IonsNamelist}) = "IONS"
QuantumESPRESSOBase.titleof(::Type{<:CellNamelist}) = "CELL"

"""
    bravais_lattice(nml::SystemNamelist)

Return a 3x3 matrix representing the Bravais lattice from `nml`.
"""
QuantumESPRESSOBase.bravais_lattice(nml::SystemNamelist) =
    bravais_lattice(nml.ibrav, nml.celldm)

function QuantumESPRESSOBase.cell_volume(nml::SystemNamelist)
    iszero(nml.ibrav) && error("`ibrav` must be non-zero to calculate the cell volume!")
    return det(bravais_lattice(nml))
end # function QuantumESPRESSOBase.cell_volume

function Setters.batchset(::VerbositySetter{:high}, template::ControlNamelist)
    return setproperties(
        template,
        verbosity = "high",
        wf_collect = true,
        tstress = true,
        tprnfor = true,
        saverho = true,
        disk_io = "high",
    )
end # function Setters.batchset
function Setters.batchset(::VerbositySetter{:low}, template::ControlNamelist)
    return setproperties(
        template,
        verbosity = "low",
        wf_collect = false,
        tstress = false,
        tprnfor = false,
        saverho = false,
        disk_io = "default",
    )
end # function Setters.batchset

end
