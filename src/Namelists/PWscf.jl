"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using LinearAlgebra: det

using Parameters: @with_kw

using QuantumESPRESSOBase
using ..Namelists: Namelist

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist

const Maybe{T} = Union{T,Nothing}

# The default values are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90.
@with_kw struct ControlNamelist <: Namelist
    calculation::String = "scf"
    title::String = " "
    verbosity::String = "low"
    restart_mode::String = "from_scratch"
    wf_collect::Bool = true
    nstep::Int = 50
    iprint::Int = 100000
    tstress::Bool = false
    tprnfor::Bool = false
    dt::Float64 = 20.0
    outdir::String = "./"
    wfcdir::String = "./"
    prefix::String = "pwscf"
    lkpoint_dir::Bool = true
    max_seconds::Float64 = 10000000.0
    etot_conv_thr::Float64 = 1e-4
    forc_conv_thr::Float64 = 1e-3
    disk_io::String = ifelse(calculation == "scf", "low", "medium")
    pseudo_dir::String = raw"$HOME/espresso/pseudo/"
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
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1282-L1369.
    @assert(calculation ∈ ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    @assert(verbosity ∈ ("high", "low", "debug", "medium", "default", "minimal"))
    @assert(restart_mode ∈ ("from_scratch", "restart"))
    @assert(nstep >= 0, "`nstep` $nstep out of range!")
    @assert(iprint >= 1, "`iprint` $iprint out of range!")
    @assert(disk_io ∈ ("high", "medium", "low", "none", "default"))
    @assert(dt >= 0, "`dt` $dt out of range!")
    @assert(max_seconds >= 0, "`max_seconds` $max_seconds out of range!")
    @assert(etot_conv_thr >= 0, "`etot_conv_thr` $etot_conv_thr out of range!")
    @assert(forc_conv_thr >= 0, "`forc_conv_thr` $forc_conv_thr out of range!")
    @assert(
        !all((gate, tefield, !dipfield)),
        "`gate` cannot be used with `tefield` if dipole correction is not active!"
    )
    @assert(
        !all((gate, dipfield, !tefield)),
        "Dipole correction is not active if `tefield = false`!"
    )
end # struct ControlNamelist

@with_kw struct SystemNamelist <: Namelist
    ibrav::Int = 1  # The default value in QE's source code is -1
    celldm::Vector{Maybe{Float64}} = [nothing]  # Must specify
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
    starting_charge::Vector{Maybe{Float64}} = []
    tot_magnetization::Float64 = -1.0
    starting_magnetization::Vector{Maybe{Float64}} = []
    ecutwfc::Float64 = 0.0
    ecutrho::Float64 = 0.0
    ecutfock::Float64 = 0.0
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
    localization_thr::Float64 = 0.0  # This is only for QE 6.4
    lda_plus_u::Bool = false
    lda_plus_u_kind::Int = 0
    Hubbard_U::Vector{Maybe{Float64}} = []
    Hubbard_J0::Vector{Maybe{Float64}} = []
    Hubbard_alpha::Vector{Maybe{Float64}} = []
    Hubbard_beta::Vector{Maybe{Float64}} = []
    # Hubbard_J::Vector{Vector{Maybe{Float64}}} = [zeros(ntyp)]  # The default value in QE's source code is just one 0.0
    starting_ns_eigenvalue::Float64 = -1.0  # It's actually a multidimensional array.
    U_projection_type::String = "atomic"
    edir::Int = 1
    emaxpos::Float64 = 0.5
    eopreg::Float64 = 0.1
    eamp::Float64 = 0.001  # The default value in QE's source code is 0.0
    angle1::Vector{Maybe{Float64}} = []
    angle2::Vector{Maybe{Float64}} = []
    constrained_magnetization::String = "none"
    fixed_magnetization::Vector{Maybe{Float64}} = zeros(3)  # The default value in QE's source code is just one 0.0
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
    london_c6::Vector{Maybe{Float64}} = []
    london_rvdw::Vector{Maybe{Float64}} = []
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
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1378-L1499.
    @assert(ibrav ∈ union(0:1:14, (-3, -5, -9, -12)))
    @assert(if ibrav != 0
        celldm[1] != 0 || A != 0  # Cannot use `iszero` to compare now!
    end, "Invalid lattice parameters (`celldm` $celldm or `A` $A)!")
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
        "`celldm` has length between 1 to 6! See `ibrav`'s doc!"
    )
    @assert(nat >= 0, "`nat` $nat is less than zero!")
    @assert(0 <= ntyp <= 10, "`ntyp` $ntyp is either less than zero or too large!")
    @assert(ntyp <= nat, "`ntyp` cannot be larger than `nat`!")
    @assert(nspin ∈ (1, 2, 4), "`nspin` $nspin out of range!")
    @assert(ecutwfc >= 0, "`ecutwfc` $ecutwfc out of range!")
    @assert(ecutrho >= 0, "`ecutrho` $ecutrho out of range!")
    @assert(ecfixed >= 0, "`ecfixed` $ecfixed out of range!")
    @assert(qcutz >= 0, "`qcutz` $qcutz out of range!")
    @assert(q2sigma >= 0, "`q2sigma` $q2sigma out of range!")
    @assert(length(starting_charge) <= ntyp)
    @assert(length(starting_magnetization) <= ntyp)
    @assert(length(Hubbard_U) <= ntyp)
    @assert(length(Hubbard_J0) <= ntyp)
    @assert(length(Hubbard_alpha) <= ntyp)
    @assert(length(Hubbard_beta) <= ntyp)
    # @assert(all(length(x) <= ntyp for x in Hubbard_J))
    @assert(length(angle1) <= ntyp)
    @assert(length(angle2) <= ntyp)
    @assert(length(fixed_magnetization) <= 3)
    @assert(length(london_c6) <= ntyp)
    @assert(length(london_rvdw) <= ntyp)
    @assert(
        exxdiv_treatment ∈ (
            "gygi-baldereschi",
            "gygi-bald",
            "g-b",
            "vcut_ws",
            "vcut_spherical",
            "none",
        ),
        "Invalid `exxdiv_treatment` $(exxdiv_treatment)!"
    )
    @assert(
        !(x_gamma_extrapolation && exxdiv_treatment ∈ ("vcut_ws", "vcut_spherical")),
        "`x_gamma_extrapolation` cannot be used with `vcut`!"
    )
end # struct SystemNamelist

@with_kw struct ElectronsNamelist <: Namelist
    electron_maxstep::Int = 100
    scf_must_converge::Bool = true
    conv_thr::Float64 = 1e-6
    adaptive_thr::Bool = false
    conv_thr_init::Float64 = 1e-3
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
    efield_cart::Vector{Maybe{Float64}} = zeros(3)
    efield_phase::String = "none"
    startingpot::String = "atomic"  # This depends on `calculation`
    startingwfc::String = "atomic+random"  # This depends on `calculation`
    tqr::Bool = false
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1508-L1543.
    @assert(
        mixing_mode ∈ ("plain", "TF", "local-TF"),
        "Invalid `mixing_mode` $(mixing_mode)!"
    )
    @assert(
        diagonalization ∈ ("david", "cg", "cg-serial", "david-serial"),
        "Invalid `diagonalization` $(diagonalization)!"
    )
    @assert(
        efield_phase ∈ ("read", "write", "none"),
        "Invalid `efield_phase` $(efield_phase)!"
    )
    @assert(startingpot ∈ ("atomic", "file"), "Invalid `startingpot` $(startingpot)!")
    @assert(
        startingwfc ∈ ("atomic", "atomic+random", "random", "file"),
        "Invalid `startingwfc` $(startingwfc)!"
    )
end # struct ElectronsNamelist

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
    trust_radius_min::Float64 = 1e-3  # The default value in QE's source code is 0.0001
    trust_radius_ini::Float64 = 0.5
    w_1::Float64 = 0.01
    w_2::Float64 = 0.5
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1552-L1585.
    @assert(
        ion_dynamics ∈ (
            "none",
            "bfgs",
            "damp",
            "verlet",
            "langevin",
            "langevin-smc",
            "beeman",
        ),
        "Invalid `ion_dynamics` $(ion_dynamics)!"
    )
    @assert(
        ion_positions ∈ ("default", "from_input"),
        "Invalid `ion_position` $(ion_positions)!"
    )
    @assert(
        pot_extrapolation ∈ ("none", "atomic", "first_order", "second_order"),
        "Invalid `pot_extrapolation` $(pot_extrapolation)!"
    )
    @assert(
        wfc_extrapolation ∈ ("none", "first_order", "second_order"),
        "Invalid `wfc_extrapolation` $(wfc_extrapolation)!"
    )
    @assert(
        ion_temperature ∈ (
            "rescaling",
            "rescale-v",
            "rescale-T",
            "reduce-T",
            "berendsen",
            "andersen",
            "initial",
            "not_controlled",
        ),
        "Invalid `ion_temperature` $(ion_temperature)!"
    )
    @assert(tempw > 0, "`tempw` $tempw out of range!")
end # struct IonsNamelist

@with_kw struct CellNamelist <: Namelist
    cell_dynamics::String = "none"
    press::Float64 = 0.0
    wmass::Float64 = 0.0
    cell_factor::Float64 = 0.0
    press_conv_thr::Float64 = 0.5
    cell_dofree::String = "all"
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1596-L1625.
    @assert(
        cell_dynamics ∈ ("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"),
        "Invalid `cell_dynamics` $(cell_dynamics)!"
    )
    @assert(wmass >= 0, "`wmass` $wmass out of range!")
    @assert(
        cell_dofree ∈ (
            "all",
            "ibrav",
            "x",
            "y",
            "z",
            "xy",
            "xz",
            "yz",
            "xyz",
            "shape",
            "volume",
            "2Dxy",
            "2Dshape",
            "epitaxial_ab",  # New in 6.4
            "epitaxial_ac",  # New in 6.4
            "epitaxial_bc",  # New in 6.4
        ),
        "Invalid `cell_dofree` $(cell_dofree)!"
    )
end # struct CellNamelist

# The following default values are picked from `<QE source>/PP/src/dos.f90`
@with_kw struct DosNamelist <: Namelist
    prefix::String = "pwscf"
    outdir::String = "./"
    ngauss::Int = 0
    degauss::Float64 = 0.0
    Emin::Float64 = -1000000.0
    Emax::Float64 = 1000000.0
    DeltaE::Float64 = 0.01
    fildos::String = "$(prefix).dos"
    @assert(ngauss ∈ (0, 1, -1, -99), "Invalid `ngauss` $(ngauss)!")
end # struct DosNamelist

# The following default values are picked from `<QE source>/PP/src/bands.f90`
@with_kw struct BandsNamelist <: Namelist
    prefix::String = "pwscf"
    outdir::String = "./"
    filband::String = "bands.out"
    spin_component::Int = 1
    lsigma::Vector{Maybe{Bool}} = falses(3)  # The default value in QE's source code is just one `false`
    lp::Bool = false
    filp::String = "p_avg.dat"
    lsym::Bool = true
    no_overlap::Bool = true
    plot_2d::Bool = false
    firstk::Int = 0
    lastk::Int = 10000000
    @assert(spin_component ∈ 1:2, "Invalid `spin_component` $(spin_component)!")
end # struct BandsNamelist

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

function QuantumESPRESSOBase.cell_volume(nml::SystemNamelist)
    iszero(nml.ibrav) && error("`ibrav` must be non-zero to calculate the cell volume!")
    return det(bravais_lattice(nml))
end # function QuantumESPRESSOBase.cell_volume

end
