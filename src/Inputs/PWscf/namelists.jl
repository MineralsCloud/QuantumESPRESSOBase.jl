# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}

# The default values are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90.
"""
    ControlNamelist <: Namelist
    ControlNamelist(calculation, title, verbosity, restart_mode, wf_collect, nstep, iprint, tstress, tprnfor, dt, outdir, wfcdir, prefix, lkpoint_dir, max_seconds, etot_conv_thr, forc_conv_thr, disk_io, pseudo_dir, tefield, dipfield, lelfield, nberrycyc, lorbm, lberry, gdir, nppstr, lfcpopt, gate)
    ControlNamelist(; kwargs...)
    ControlNamelist(::ControlNamelist, kwargs...)
    ControlNamelist(::ControlNamelist, dict)

Represent the `CONTROL` namelist of `pw.x`.
"""
struct ControlNamelist <: Namelist
    calculation::String
    title::String
    verbosity::String
    restart_mode::String
    wf_collect::Bool
    nstep::UInt
    iprint::UInt
    tstress::Bool
    tprnfor::Bool
    dt::Float64
    outdir::String
    wfcdir::String
    prefix::String
    lkpoint_dir::Bool
    max_seconds::Float64
    etot_conv_thr::Float64
    forc_conv_thr::Float64
    disk_io::String
    pseudo_dir::String
    tefield::Bool
    dipfield::Bool
    lelfield::Bool
    nberrycyc::UInt
    lorbm::Bool
    lberry::Bool
    gdir::UInt
    nppstr::UInt
    lfcpopt::Bool
    gate::Bool
end # struct ControlNamelist
function ControlNamelist(;
    calculation = "scf",
    title = " ",
    verbosity = "low",
    restart_mode = "from_scratch",
    wf_collect = true,
    nstep = 50,
    iprint = 100000,
    tstress = false,
    tprnfor = false,
    dt = 20.0,
    outdir = "./",
    wfcdir = "./",
    prefix = "pwscf",
    lkpoint_dir = true,
    max_seconds = 10000000.0,
    etot_conv_thr = 1e-4,
    forc_conv_thr = 1e-3,
    disk_io = ifelse(calculation == "scf", "low", "medium"),
    pseudo_dir = raw"$HOME/espresso/pseudo/",
    tefield = false,
    dipfield = false,
    lelfield = false,
    nberrycyc = 1,
    lorbm = false,
    lberry = false,
    gdir = 1,  # The QE default value is `0`!
    nppstr = 0,
    lfcpopt = false,
    gate = false,
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1282-L1369.
    @argcheck calculation in ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
    @argcheck verbosity in ("high", "low", "debug", "medium", "default", "minimal")
    @argcheck restart_mode in ("from_scratch", "restart")
    @argcheck iprint >= 1
    @argcheck disk_io in ("high", "medium", "low", "none", "default")
    @argcheck dt >= 0
    # @argcheck !(lkpoint_dir && wf_collect) "`lkpoint_dir` currently doesn't work together with `wf_collect`!"
    @argcheck max_seconds >= 0
    @argcheck etot_conv_thr >= 0
    @argcheck forc_conv_thr >= 0
    @argcheck gdir in 1:3
    @argcheck !all((gate, tefield, !dipfield)) "`gate` cannot be used with `tefield` if dipole correction is not active!"
    @argcheck !all((gate, dipfield, !tefield)) "dipole correction is not active if `tefield = false`!"
    return ControlNamelist(
        calculation,
        title,
        verbosity,
        restart_mode,
        wf_collect,
        nstep,
        iprint,
        tstress,
        tprnfor,
        dt,
        outdir,
        wfcdir,
        prefix,
        lkpoint_dir,
        max_seconds,
        etot_conv_thr,
        forc_conv_thr,
        disk_io,
        pseudo_dir,
        tefield,
        dipfield,
        lelfield,
        nberrycyc,
        lorbm,
        lberry,
        gdir,
        nppstr,
        lfcpopt,
        gate,
    )
end
ControlNamelist(nml::ControlNamelist; kwargs...) = setproperties(nml, kwargs...)
ControlNamelist(nml::ControlNamelist, t::NamedTuple) = setproperties(nml, t)
ControlNamelist(nml::ControlNamelist, dict::AbstractDict) = setproperties(nml, dict)

xmldir(nml::ControlNamelist) = expanduser(joinpath(nml.outdir, nml.prefix * ".save"))
wfcfiles(nml::ControlNamelist, n = 1) =
    [joinpath(xmldir(nml), nml.prefix * ".wfc$i") for i = 1:n]

"""
    SystemNamelist <: Namelist
    SystemNamelist(ibrav, celldm, A, B, C, cosAB, cosAC, cosBC, nat, ntyp, nbnd, tot_charge, starting_charge, tot_magnetization, starting_magnetization, ecutwfc, ecutrho, ecutfock, nr1, nr2, nr3, nr1s, nr2s, nr3s, nosym, nosym_evc, noinv, no_t_rev, force_symmorphic, use_all_frac, occupations, one_atom_occupations, starting_spin_angle, degauss, smearing, nspin, noncolin, ecfixed, qcutz, q2sigma, input_dft, exx_fraction, screening_parameter, exxdiv_treatment, x_gamma_extrapolation, ecutvcut, nqx1, nqx2, nqx3, localization_thr, lda_plus_u, lda_plus_u_kind, Hubbard_U, Hubbard_J0, Hubbard_alpha, Hubbard_beta, starting_ns_eigenvalue, U_projection_type, edir, emaxpos, eopreg, eamp, angle1, angle2, constrained_magnetization, fixed_magnetization, lambda, report, lspinorb, assume_isolated, esm_bc, esm_w, esm_efield, esm_nfit, fcp_mu, vdw_corr, london, london_s6, london_c6, london_rvdw, london_rcut, ts_vdw_econv_thr, ts_vdw_isolated, xdm, xdm_a1, xdm_a2, space_group, uniqueb, origin_choice, rhombohedral, zgate, relaxz, block, block_1, block_2, block_height)
    SystemNamelist(; kwargs...)
    SystemNamelist(::SystemNamelist, kwargs...)
    SystemNamelist(::SystemNamelist, dict)

Represent the `SYSTEM` namelist of `pw.x`.
"""
struct SystemNamelist <: Namelist
    ibrav::UInt
    celldm::Vector{Maybe{Float64}}
    A::Float64
    B::Float64
    C::Float64
    cosAB::Float64
    cosAC::Float64
    cosBC::Float64
    nat::UInt
    ntyp::UInt
    nbnd::UInt
    tot_charge::Float64
    starting_charge::Vector{Maybe{Float64}}
    tot_magnetization::Float64
    starting_magnetization::Vector{Maybe{Float64}}
    ecutwfc::Float64
    ecutrho::Float64
    ecutfock::Float64
    nr1::UInt
    nr2::UInt
    nr3::UInt
    nr1s::UInt
    nr2s::UInt
    nr3s::UInt
    nosym::Bool
    nosym_evc::Bool
    noinv::Bool
    no_t_rev::Bool
    force_symmorphic::Bool
    use_all_frac::Bool
    occupations::String
    one_atom_occupations::Bool
    starting_spin_angle::Bool
    degauss::Float64
    smearing::String
    nspin::UInt
    noncolin::Bool
    ecfixed::Float64
    qcutz::Float64
    q2sigma::Float64
    input_dft::String
    exx_fraction::Float64
    screening_parameter::Float64
    exxdiv_treatment::String
    x_gamma_extrapolation::Bool
    ecutvcut::Float64
    nqx1::UInt
    nqx2::UInt
    nqx3::UInt
    localization_thr::Float64
    lda_plus_u::Bool
    lda_plus_u_kind::UInt
    Hubbard_U::Vector{Maybe{Float64}}
    Hubbard_J0::Vector{Maybe{Float64}}
    Hubbard_alpha::Vector{Maybe{Float64}}
    Hubbard_beta::Vector{Maybe{Float64}}
    # Hubbard_J::Vector{Vector{Maybe{Float64}}}
    starting_ns_eigenvalue::Float64
    U_projection_type::String
    edir::UInt
    emaxpos::Float64
    eopreg::Float64
    eamp::Float64
    angle1::Vector{Maybe{Float64}}
    angle2::Vector{Maybe{Float64}}
    constrained_magnetization::String
    fixed_magnetization::Vector{Maybe{Float64}}
    lambda::Float64
    report::UInt
    lspinorb::Bool
    assume_isolated::String
    esm_bc::String
    esm_w::Float64
    esm_efield::Float64
    esm_nfit::UInt
    fcp_mu::Float64
    vdw_corr::String
    london::Bool
    london_s6::Float64
    london_c6::Vector{Maybe{Float64}}
    london_rvdw::Vector{Maybe{Float64}}
    london_rcut::Float64
    ts_vdw_econv_thr::Float64
    ts_vdw_isolated::Bool
    xdm::Bool
    xdm_a1::Float64
    xdm_a2::Float64
    space_group::UInt
    uniqueb::Bool
    origin_choice::UInt
    rhombohedral::Bool
    zgate::Float64
    relaxz::Bool
    block::Bool
    block_1::Float64
    block_2::Float64
    block_height::Float64
end # struct SystemNamelist
function SystemNamelist(;
    ibrav,
    celldm = zeros(6),  # Must specify
    A = 0.0,
    B = 0.0,
    C = 0.0,
    cosAB = 0.0,
    cosAC = 0.0,
    cosBC = 0.0,
    nat = 0,
    ntyp = 0,
    nbnd = 0,
    tot_charge = 0.0,
    starting_charge = [],
    tot_magnetization = -1.0,
    starting_magnetization = [],
    ecutwfc = 0.0,
    ecutrho = 0.0,
    ecutfock = 0.0,
    nr1 = 0,
    nr2 = 0,
    nr3 = 0,
    nr1s = 0,
    nr2s = 0,
    nr3s = 0,
    nosym = false,
    nosym_evc = false,
    noinv = false,
    no_t_rev = false,
    force_symmorphic = false,
    use_all_frac = false,
    occupations = "fixed",
    one_atom_occupations = false,
    starting_spin_angle = false,
    degauss = 0.0,
    smearing = "gaussian",
    nspin = 1,
    noncolin = false,
    ecfixed = 0.0,
    qcutz = 0.0,
    q2sigma = 0.1,  # The default value in QE's source code is 0.01
    input_dft = "none",
    exx_fraction = 0.25,
    screening_parameter = 0.106,
    exxdiv_treatment = "gygi-baldereschi",
    x_gamma_extrapolation = true,
    ecutvcut = 0.0,
    nqx1 = 1,
    nqx2 = 1,
    nqx3 = 1,
    localization_thr = 0.0,  # This is only for QE 6.4
    lda_plus_u = false,
    lda_plus_u_kind = 0,
    Hubbard_U = [],
    Hubbard_J0 = [],
    Hubbard_alpha = [],
    Hubbard_beta = [],
    # Hubbard_J = [zeros(ntyp)]  ,  # The default value in QE's source code is just one 0.0
    starting_ns_eigenvalue = -1.0,  # It's actually a multidimensional array.
    U_projection_type = "atomic",
    edir = 1,
    emaxpos = 0.5,
    eopreg = 0.1,
    eamp = 0.001,  # The default value in QE's source code is 0.0
    angle1 = [],
    angle2 = [],
    constrained_magnetization = "none",
    fixed_magnetization = zeros(3),  # The default value in QE's source code is just one 0.0
    lambda = 1.0,
    report = 100,
    lspinorb = false,
    assume_isolated = "none",
    esm_bc = "pbc",
    esm_w = 0.0,
    esm_efield = 0.0,
    esm_nfit = 4,
    fcp_mu = 0.0,
    vdw_corr = "none",
    london = false,
    london_s6 = 0.75,
    london_c6 = [],
    london_rvdw = [],
    london_rcut = 200.0,
    ts_vdw_econv_thr = 1e-06,
    ts_vdw_isolated = false,
    xdm = false,
    xdm_a1 = 0.6836,  # The default value in QE's source code is 0.0
    xdm_a2 = 1.5045,  # The default value in QE's source code is 0.0
    space_group = 0,
    uniqueb = false,
    origin_choice = 1,
    rhombohedral = true,
    zgate = 0.5,
    relaxz = false,
    block = false,
    block_1 = 0.45,
    block_2 = 0.55,
    block_height = 0.1,  # The default value in QE's source code is 0.0
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1378-L1499.
    @argcheck ibrav in union(0:1:14, (-3, -5, -9, 91, -12, -13))
    @argcheck ntyp <= 10 "`ntyp` $ntyp is larger than 10!"
    @argcheck ntyp <= nat
    @argcheck smearing in (
        "gaussian",
        "gauss",
        "methfessel-paxton",
        "m-p",
        "mp",
        "marzari-vanderbilt",
        "cold",
        "m-v",
        "mv",
        "fermi-dirac",
        "f-d",
        "fd",
    )
    @argcheck nspin in (1, 2, 4)
    @argcheck ecutwfc >= 0
    @argcheck ecutrho >= 0
    @argcheck ecfixed >= 0
    @argcheck qcutz >= 0
    @argcheck q2sigma >= 0
    @argcheck lda_plus_u_kind in 0:1
    @argcheck edir in 1:3
    @argcheck origin_choice in 1:2
    @argcheck length(starting_charge) <= ntyp
    @argcheck length(starting_magnetization) <= ntyp
    @argcheck length(Hubbard_U) <= ntyp
    @argcheck length(Hubbard_J0) <= ntyp
    @argcheck length(Hubbard_alpha) <= ntyp
    @argcheck length(Hubbard_beta) <= ntyp
    # @argcheck all(length(x) <= ntyp for x in Hubbard_J)
    @argcheck length(angle1) <= ntyp
    @argcheck length(angle2) <= ntyp
    @argcheck length(fixed_magnetization) <= 3
    @argcheck length(london_c6) <= ntyp
    @argcheck length(london_rvdw) <= ntyp
    @argcheck exxdiv_treatment in
              ("gygi-baldereschi", "gygi-bald", "g-b", "vcut_ws", "vcut_spherical", "none")
    @argcheck !(x_gamma_extrapolation && exxdiv_treatment in ("vcut_ws", "vcut_spherical")) "`x_gamma_extrapolation` cannot be used with `vcut`!"
    return SystemNamelist(
        ibrav,
        celldm,
        A,
        B,
        C,
        cosAB,
        cosAC,
        cosBC,
        nat,
        ntyp,
        nbnd,
        tot_charge,
        starting_charge,
        tot_magnetization,
        starting_magnetization,
        ecutwfc,
        ecutrho,
        ecutfock,
        nr1,
        nr2,
        nr3,
        nr1s,
        nr2s,
        nr3s,
        nosym,
        nosym_evc,
        noinv,
        no_t_rev,
        force_symmorphic,
        use_all_frac,
        occupations,
        one_atom_occupations,
        starting_spin_angle,
        degauss,
        smearing,
        nspin,
        noncolin,
        ecfixed,
        qcutz,
        q2sigma,
        input_dft,
        exx_fraction,
        screening_parameter,
        exxdiv_treatment,
        x_gamma_extrapolation,
        ecutvcut,
        nqx1,
        nqx2,
        nqx3,
        localization_thr,
        lda_plus_u,
        lda_plus_u_kind,
        Hubbard_U,
        Hubbard_J0,
        Hubbard_alpha,
        Hubbard_beta,
        starting_ns_eigenvalue,
        U_projection_type,
        edir,
        emaxpos,
        eopreg,
        eamp,
        angle1,
        angle2,
        constrained_magnetization,
        fixed_magnetization,
        lambda,
        report,
        lspinorb,
        assume_isolated,
        esm_bc,
        esm_w,
        esm_efield,
        esm_nfit,
        fcp_mu,
        vdw_corr,
        london,
        london_s6,
        london_c6,
        london_rvdw,
        london_rcut,
        ts_vdw_econv_thr,
        ts_vdw_isolated,
        xdm,
        xdm_a1,
        xdm_a2,
        space_group,
        uniqueb,
        origin_choice,
        rhombohedral,
        zgate,
        relaxz,
        block,
        block_1,
        block_2,
        block_height,
    )
end
SystemNamelist(nml::SystemNamelist; kwargs...) = setproperties(nml, kwargs...)
SystemNamelist(nml::SystemNamelist, t::NamedTuple) = setproperties(nml, t)
SystemNamelist(nml::SystemNamelist, dict::AbstractDict) = setproperties(nml, dict)

"""
    ElectronsNamelist <: Namelist
    ElectronsNamelist(electron_maxstep, scf_must_converge, conv_thr, adaptive_thr, conv_thr_init, conv_thr_multi, mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns, diagonalization, ortho_para, diago_thr_init, diago_cg_maxiter, diago_david_ndim, diago_full_acc, efield, efield_cart, efield_phase, startingpot, startingwfc, tqr)
    ElectronsNamelist(; kwargs...)
    ElectronsNamelist(::ElectronsNamelist, kwargs...)
    ElectronsNamelist(::ElectronsNamelist, dict)

Represent the `ELECTRONS` namelist of `pw.x`.
"""
struct ElectronsNamelist <: Namelist
    electron_maxstep::UInt
    scf_must_converge::Bool
    conv_thr::Float64
    adaptive_thr::Bool
    conv_thr_init::Float64
    conv_thr_multi::Float64
    mixing_mode::String
    mixing_beta::Float64
    mixing_ndim::UInt
    mixing_fixed_ns::UInt
    diagonalization::String
    ortho_para::UInt
    diago_thr_init::Float64
    diago_cg_maxiter::UInt
    diago_david_ndim::UInt
    diago_full_acc::Bool
    efield::Float64
    efield_cart::Vector{Maybe{Float64}}
    efield_phase::String
    startingpot::String  # This depends on `calculation`
    startingwfc::String  # This depends on `calculation`
    tqr::Bool
end # struct ElectronsNamelist
function ElectronsNamelist(;
    electron_maxstep = 100,
    scf_must_converge = true,
    conv_thr = 1e-6,
    adaptive_thr = false,
    conv_thr_init = 1e-3,
    conv_thr_multi = 0.1,
    mixing_mode = "plain",
    mixing_beta = 0.7,
    mixing_ndim = 8,
    mixing_fixed_ns = 0,
    diagonalization = "david",
    ortho_para = 0,
    diago_thr_init = 0.0,
    diago_cg_maxiter = 20,
    diago_david_ndim = 4,
    diago_full_acc = false,
    efield = 0.0,
    efield_cart = zeros(3),
    efield_phase = "none",
    startingpot = "atomic",  # This depends on `calculation`
    startingwfc = "atomic+random",  # This depends on `calculation`
    tqr = false,
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1508-L1543.
    @argcheck mixing_mode in ("plain", "TF", "local-TF")
    @argcheck diagonalization in ("david", "cg", "cg-serial", "david-serial", "ppcg")  # Different from docs
    @argcheck efield_phase in ("read", "write", "none")
    @argcheck startingpot in ("atomic", "file")
    @argcheck startingwfc in ("atomic", "atomic+random", "random", "file")
    return ElectronsNamelist(
        electron_maxstep,
        scf_must_converge,
        conv_thr,
        adaptive_thr,
        conv_thr_init,
        conv_thr_multi,
        mixing_mode,
        mixing_beta,
        mixing_ndim,
        mixing_fixed_ns,
        diagonalization,
        ortho_para,
        diago_thr_init,
        diago_cg_maxiter,
        diago_david_ndim,
        diago_full_acc,
        efield,
        efield_cart,
        efield_phase,
        startingpot,
        startingwfc,
        tqr,
    )
end
ElectronsNamelist(nml::ElectronsNamelist; kwargs...) = setproperties(nml, kwargs...)
ElectronsNamelist(nml::ElectronsNamelist, t::NamedTuple) = setproperties(nml, t)
ElectronsNamelist(nml::ElectronsNamelist, dict::AbstractDict) = setproperties(nml, dict)

"""
    IonsNamelist <: Namelist
    IonsNamelist(ion_dynamics, ion_positions, pot_extrapolation, wfc_extrapolation, remove_rigid_rot, ion_temperature, tempw, tolp, delta_t, nraise, refold_pos, upscale, bfgs_ndim, trust_radius_max, trust_radius_min, trust_radius_ini, w_1, w_2)
    IonsNamelist(; kwargs...)
    IonsNamelist(::IonsNamelist, kwargs...)
    IonsNamelist(::IonsNamelist, dict)

Represent the `IONS` namelist of `pw.x`.

Input this namelist only if `calculation` is `"relax"`, `"md"`, `"vc-relax"`, or `"vc-md"`.
"""
struct IonsNamelist <: Namelist
    ion_dynamics::String
    ion_positions::String
    pot_extrapolation::String
    wfc_extrapolation::String
    remove_rigid_rot::Bool
    ion_temperature::String
    tempw::Float64
    tolp::Float64
    delta_t::Float64
    nraise::UInt
    refold_pos::Bool
    upscale::Float64
    bfgs_ndim::UInt
    trust_radius_max::Float64
    trust_radius_min::Float64  # The default value in QE's source code is 0.0001
    trust_radius_ini::Float64
    w_1::Float64
    w_2::Float64
end # struct IonsNamelist
function IonsNamelist(;
    ion_dynamics = "none",
    ion_positions = "default",
    pot_extrapolation = "atomic",
    wfc_extrapolation = "none",
    remove_rigid_rot = false,
    ion_temperature = "not_controlled",
    tempw = 300.0,
    tolp = 100.0,
    delta_t = 1.0,
    nraise = 1,
    refold_pos = false,
    upscale = 100.0,
    bfgs_ndim = 1,
    trust_radius_max = 0.8,
    trust_radius_min = 1e-3,  # The default value in QE's source code is 0.0001
    trust_radius_ini = 0.5,
    w_1 = 0.01,
    w_2 = 0.5,
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1552-L1585.
    @argcheck ion_dynamics in
              ("none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman")
    @argcheck ion_positions in ("default", "from_input")
    @argcheck pot_extrapolation in ("none", "atomic", "first_order", "second_order")
    @argcheck wfc_extrapolation in ("none", "first_order", "second_order")
    @argcheck ion_temperature in (
        "rescaling",
        "rescale-v",
        "rescale-T",
        "reduce-T",
        "berendsen",
        "andersen",
        "initial",
        "not_controlled",
    )
    @argcheck tempw > 0
    return IonsNamelist(
        ion_dynamics,
        ion_positions,
        pot_extrapolation,
        wfc_extrapolation,
        remove_rigid_rot,
        ion_temperature,
        tempw,
        tolp,
        delta_t,
        nraise,
        refold_pos,
        upscale,
        bfgs_ndim,
        trust_radius_max,
        trust_radius_min,
        trust_radius_ini,
        w_1,
        w_2,
    )
end
IonsNamelist(nml::IonsNamelist; kwargs...) = setproperties(nml, kwargs...)
IonsNamelist(nml::IonsNamelist, t::NamedTuple) = setproperties(nml, t)
IonsNamelist(nml::IonsNamelist, dict::AbstractDict) = setproperties(nml, dict)

"""
    CellNamelist <: Namelist
    CellNamelist(cell_dynamics, press, wmass, cell_factor, press_conv_thr, cell_dofree)
    CellNamelist(; kwargs...)
    CellNamelist(::CellNamelist, kwargs...)
    CellNamelist(::CellNamelist, dict)

Represent the `CELL` namelist of `pw.x`.

Input this namelist only if `calculation` is `"vc-relax"` or `"vc-md"`.
"""
struct CellNamelist <: Namelist
    cell_dynamics::String
    press::Float64
    wmass::Float64
    cell_factor::Float64
    press_conv_thr::Float64
    cell_dofree::String
end # struct CellNamelist
function CellNamelist(;
    cell_dynamics = "none",
    press = 0.0,
    wmass = 0.0,
    cell_factor = 0.0,
    press_conv_thr = 0.5,
    cell_dofree = "all",
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1596-L1625.
    @argcheck cell_dynamics in ("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w")
    @argcheck wmass >= 0
    @argcheck cell_dofree in (
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
    )
    return CellNamelist(
        cell_dynamics,
        press,
        wmass,
        cell_factor,
        press_conv_thr,
        cell_dofree,
    )
end
CellNamelist(nml::CellNamelist; kwargs...) = setproperties(nml, kwargs...)
CellNamelist(nml::CellNamelist, t::NamedTuple) = setproperties(nml, t)
CellNamelist(nml::CellNamelist, dict::AbstractDict) = setproperties(nml, dict)

# The following default values are picked from `<QE source>/PP/src/dos.f90`
"""
    DosNamelist <: Namelist

Represent the `DOS` namelist of `dos.x`.
"""
struct DosNamelist <: Namelist
    prefix::String
    outdir::String
    ngauss::Int
    degauss::Float64
    Emin::Float64
    Emax::Float64
    DeltaE::Float64
    fildos::String
end # struct DosNamelist
function DosNamelist(;
    prefix = "pwscf",
    outdir = "./",
    ngauss = 0,
    degauss = 0.0,
    Emin = -1000000.0,
    Emax = 1000000.0,
    DeltaE = 0.01,
    fildos = "$(prefix).dos",
)
    @argcheck ngauss in (0, 1, -1, -99)
    return DosNamelist(prefix, outdir, ngauss, degauss, Emin, Emax, DeltaE, fildos)
end
DosNamelist(nml::DosNamelist; kwargs...) = setproperties(nml, kwargs...)
DosNamelist(nml::DosNamelist, t::NamedTuple) = setproperties(nml, t)
DosNamelist(nml::DosNamelist, dict::AbstractDict) = setproperties(nml, dict)

# The following default values are picked from `<QE source>/PP/src/bands.f90`
"""
    BandsNamelist <: Namelist

Represent the `BANDS` namelist of `bands.x`.
"""
struct BandsNamelist <: Namelist
    prefix::String
    outdir::String
    filband::String
    spin_component::Int
    lsigma::Vector{Maybe{Bool}}  # The default value in QE's source code is just one `false`
    lp::Bool
    filp::String
    lsym::Bool
    no_overlap::Bool
    plot_2d::Bool
    firstk::UInt
    lastk::UInt
end # struct BandsNamelist
function BandsNamelist(;
    prefix = "pwscf",
    outdir = "./",
    filband = "bands.out",
    spin_component = 1,
    lsigma = falses(3),  # The default value in QE's source code is just one `false`
    lp = false,
    filp = "p_avg.dat",
    lsym = true,
    no_overlap = true,
    plot_2d = false,
    firstk = 0,
    lastk = 10000000,
)
    @argcheck spin_component in 1:2
    return BandsNamelist(
        prefix,
        outdir,
        filband,
        spin_component,
        lsigma,
        lp,
        filp,
        lsym,
        no_overlap,
        plot_2d,
        firstk,
        lastk,
    )
end
BandsNamelist(nml::BandsNamelist; kwargs...) = setproperties(nml, kwargs...)
BandsNamelist(nml::BandsNamelist, t::NamedTuple) = setproperties(nml, t)
BandsNamelist(nml::BandsNamelist, dict::AbstractDict) = setproperties(nml, dict)

"""
    set_verbosity(template::ControlNamelist, verbosity)

Return a modified `ControlNamelist`, with verbosity set.
"""
function set_verbosity(control::ControlNamelist, verbosity)
    if verbosity == "high"
        @set! control.verbosity = "high"
        @set! control.wf_collect = true
        @set! control.tstress = true
        @set! control.tprnfor = true
        @set! control.disk_io = "high"
    elseif verbosity == "low"
        @set! control.verbosity = "low"
        @set! control.wf_collect = false
        @set! control.tstress = false
        @set! control.tprnfor = false
        @set! control.disk_io = "low"
    else
        error("unknown `verbosity` `$verbosity` specified!")
    end
    return control
end # function set_verbosity

"""
    set_temperature(system::SystemNamelist, temperature)

Return a modified `SystemNamelist`, with finite temperature set.

!!! warning
    Can be used with(out) units. If no unit is given, "Ry" is chosen.
"""
function set_temperature(system::SystemNamelist, temperature)
    @set! system.occupations = "smearing"
    @set! system.smearing = "fermi-dirac"
    @set! system.degauss = _set_temperature(temperature)
    return system
end # function set_temperature
function _set_temperature(temperature::AbstractQuantity)
    u = upreferred(unit(temperature))
    if u == u"Ry"
        return temperature
    elseif u == u"kg*m^2*s^-2"  # u"hartree", u"J", u"eV", etc..
        return ustrip(u"Ry", temperature)
    elseif u == u"s^-1"  # u"Hz", u"THz", ...
        return ustrip(u"Hz", temperature) / 6579683920502000.0 * 2
    elseif u == u"K"  # u"K", u"mK", u"μK", ...
        return ustrip(u"K", temperature) / 315775.02480407 * 2
    elseif u == u"m^-1"  # u"m^-1", u"cm^-1", ...
        return ustrip(u"m^-1", temperature) / 21947463.13632 * 2
    elseif u == u"kg"  # u"kg", u"g", ...
        return ustrip(u"kg", temperature) / 4.8508702095432e-35 * 2
    else
        error("unknown unit given!")
    end
end # function _set_temperature
_set_temperature(temperature::Real) = temperature  # Ry
