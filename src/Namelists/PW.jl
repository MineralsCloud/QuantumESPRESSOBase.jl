"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists: Namelist

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist

# The following default values are picked from `<QE source>/Modules/read_namelists.f90`
@with_kw struct ControlNamelist <: Namelist
    calculation::String = "scf"
    title::String = " "
    verbosity::String = "low"
    restart_mode::String = "from_scratch"
    wf_collect::Bool = true
    nstep::Int = 1
    iprint::Int = 1
    tstress::Bool = false
    tprnfor::Bool = false
    dt::Float64 = 20.0
    outdir::String = "./"
    wfcdir::String = "./"
    prefix::String = "pwscf"
    lkpoint_dir::Bool = true
    max_seconds::Float64 = 10000000.0
    etot_conv_thr::Float64 = 0.0001
    forc_conv_thr::Float64 = 0.001
    disk_io::String = "medium"
    pseudo_dir::String = raw"$HOME/espresso/pseudo/"
    tefield::Bool = false
    dipfield::Bool = false
    lelfield::Bool = false
    nberrycyc::Int = 1
    lorbm::Bool = false
    lberry::Bool = false
    gdir::Int = 1
    nppstr::Int = 1
    lfcpopt::Bool = false
    gate::Bool = false
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
    starting_charge::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    tot_magnetization::Float64 = -1.0
    starting_magnetization::Vector{Union{Missing, Float64}} = ones(ntyp)
    ecutwfc::Float64 = 0.0
    ecutrho::Float64 = 0.0
    ecutfock::Float64 = -1.0
    nr1::Int = 0
    nr2::Int = 0
    nr3::Int = 0
    nr1s::Int = 0
    nr2s::Int = 0
    nr3s::Int = 0
    nosym::Bool = false
    nosym_evc::Bool = false
    noinv::Bool = false
    no_t_rev::Bool = false
    force_symmorphic::Bool = false
    use_all_frac::Bool = false
    occupations::String = "fixed"
    one_atom_occupations::Bool = false
    starting_spin_angle::Bool = false
    degauss::Float64 = 0.0
    smearing::String = "gaussian"
    nspin::Int = 1
    noncolin::Bool = false
    ecfixed::Float64 = 0.0
    qcutz::Float64 = 0.0
    q2sigma::Float64 = 0.1  # The default value in QE's source code is 0.01
    input_dft::String = "none"
    exx_fraction::Float64 = 0.25
    screening_parameter::Float64 = 0.106
    exxdiv_treatment::String = "gygi-baldereschi"
    x_gamma_extrapolation::Bool = true
    ecutvcut::Float64 = 0.0
    nqx1::Int = 1
    nqx2::Int = 1
    nqx3::Int = 1
    lda_plus_u::Bool = false
    lda_plus_u_kind::Int = 0
    Hubbard_U::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    Hubbard_J0::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    Hubbard_alpha::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    Hubbard_beta::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    Hubbard_J::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    starting_ns_eigenvalue::Float64 = -1.0
    U_projection_type::String = "atomic"
    edir::Int = 1
    emaxpos::Float64 = 0.5
    eopreg::Float64 = 0.1
    eamp::Float64 = 0.001  # The default value in QE's source code is 0.0
    angle1::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    angle2::Vector{Union{Missing, Float64}} = zeros(ntyp)  # The default value in QE's source code is just one 0.0
    constrained_magnetization::String = "none"
    fixed_magnetization::Vector{Union{Missing, Float64}} = zeros(3)  # The default value in QE's source code is just one 0.0
    lambda::Float64 = 1.0
    report::Int = 100
    lspinorb::Bool = false
    assume_isolated::String = "none"
    esm_bc::String = "pbc"
    esm_w::Float64 = 0.0
    esm_efield::Float64 = 0.0
    esm_nfit::Int = 4
    fcp_mu::Float64 = 0.0
    vdw_corr::String = "none"
    london::Bool = false
    london_s6::Float64 = 0.75
    london_c6::Vector{Union{Missing, Float64}} = -ones(ntyp)  # The default value in QE's source code is just one -1.0
    london_rvdw::Vector{Union{Missing, Float64}} = -ones(ntyp)  # The default value in QE's source code is just one -1.0
    london_rcut::Float64 = 200.0
    ts_vdw_econv_thr::Float64 = 1e-06
    ts_vdw_isolated::Bool = false
    xdm::Bool = false
    xdm_a1::Float64 = 0.6836  # The default value in QE's source code is 0.0
    xdm_a2::Float64 = 1.5045  # The default value in QE's source code is 0.0
    space_group::Int = 0
    uniqueb::Bool = false
    origin_choice::Int = 1
    rhombohedral::Bool = true
    zgate::Float64 = 0.5
    relaxz::Bool = false
    block::Bool = false
    block_1::Float64 = 0.45
    block_2::Float64 = 0.55
    block_height::Float64 = 0.1  # The default value in QE's source code is 0.0
end  # struct SystemNamelist

@with_kw struct ElectronsNamelist <: Namelist
    electron_maxstep::Int = 100
    scf_must_converge::Bool = true
    conv_thr::Float64 = 1e-06
    adaptive_thr::Bool = false
    conv_thr_init::Float64 = 0.001
    conv_thr_multi::Float64 = 0.1
    mixing_mode::String = "plain"
    mixing_beta::Float64 = 0.7
    mixing_ndim::Int = 8
    mixing_fixed_ns::Int = 0
    diagonalization::String = "david"
    ortho_para::Int = 0
    diago_thr_init::Float64 = 0.0
    diago_cg_maxiter::Int = 20
    diago_david_ndim::Int = 4
    diago_full_acc::Bool = false
    efield::Float64 = 0.0
    efield_cart::Vector{Union{Missing, Float64}} = [0.0, 0.0, 0.0]
    efield_phase::String = "none"
    startingpot::String = "atomic"
    startingwfc::String = "atomic+random"
    tqr::Bool = false
end  # struct ElectronsNamelist

@with_kw struct IonsNamelist <: Namelist
    ion_dynamics::String = "none"
    ion_positions::String = "default"
    pot_extrapolation::String = "atomic"
    wfc_extrapolation::String = "none"
    remove_rigid_rot::Bool = false
    ion_temperature::String = "not_controlled"
    tempw::Float64 = 300.0
    tolp::Float64 = 100.0
    delta_t::Float64 = 1.0
    nraise::Int = 1
    refold_pos::Bool = false
    upscale::Float64 = 100.0
    bfgs_ndim::Int = 1
    trust_radius_max::Float64 = 0.8
    trust_radius_min::Float64 = 0.001  # The default value in QE's source code is 0.0001
    trust_radius_ini::Float64 = 0.5
    w_1::Float64 = 0.01
    w_2::Float64 = 0.5
end  # struct IonsNamelist

@with_kw struct CellNamelist <: Namelist
    cell_dynamics::String = "none"
    press::Float64 = 0.0
    wmass::Float64 = 0.0
    cell_factor::Float64 = 0.0
    press_conv_thr::Float64 = 0.5
    cell_dofree::String = "all"
end  # struct CellNamelist

end
