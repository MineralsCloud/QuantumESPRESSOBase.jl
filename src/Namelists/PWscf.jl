"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using LinearAlgebra: det

using Kaleido: @batchlens
using Parameters: @with_kw
using Setfield: set, @lens
using Unitful
using UnitfulAtomic

using Crystallography: BravaisLattice
using QuantumESPRESSOBase.Setters:
    VerbositySetter, FiniteTemperatureSetter, CalculationSetter, LensMaker
using QuantumESPRESSOBase.Namelists: Namelist

import Crystallography
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Setters

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist

const Maybe{T} = Union{T,Nothing}

# The default values are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90.
"""
    ControlNamelist <: Namelist

Represent the `CONTROL` namelist of `pw.x`.
"""
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
    @assert calculation ∈ ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
    @assert verbosity ∈ ("high", "low", "debug", "medium", "default", "minimal")
    @assert restart_mode ∈ ("from_scratch", "restart")
    @assert nstep >= 0
    @assert iprint >= 1
    @assert disk_io ∈ ("high", "medium", "low", "none", "default")
    @assert dt >= 0
    @assert max_seconds >= 0
    @assert etot_conv_thr >= 0
    @assert forc_conv_thr >= 0
    @assert(
        !all((gate, tefield, !dipfield)),
        "`gate` cannot be used with `tefield` if dipole correction is not active!"
    )
    @assert(
        !all((gate, dipfield, !tefield)),
        "Dipole correction is not active if `tefield = false`!"
    )
end # struct ControlNamelist

"""
    SystemNamelist <: Namelist

Represent the `SYSTEM` namelist of `pw.x`.
"""
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
    @assert ibrav ∈ union(0:1:14, (-3, -5, -9, 91, -12, -13))
    @assert(
        ibrav != 0 ? true : celldm[1] != 0 || A != 0,  # Cannot use `iszero` to compare now!
        "Invalid lattice parameters (`celldm` $celldm or `A` $A)!"
    )
    @assert(
        if ibrav == 14
            length(celldm) == 6
        elseif ibrav ∈ (5, -5, 12, 13)
            4 <= length(celldm) <= 6
        elseif ibrav ∈ (4, 6, 7, 8, 9, -9, 91, 10, 11)  # `91` is new from QE 6.4
            3 <= length(celldm) <= 6
        elseif ibrav == -13  # `-13` is new from QE 6.4
            5 <= length(celldm) <= 6
        else
            1 <= length(celldm) <= 6
        end,
        "`celldm` must have length between 1 to 6! See `ibrav`'s doc!"
    )
    @assert nat >= 0
    @assert(0 <= ntyp <= 10, "`ntyp` $ntyp is either less than zero or too large!")
    @assert ntyp <= nat
    @assert(
        smearing ∈ (
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
    )
    @assert nspin ∈ (1, 2, 4)
    @assert ecutwfc >= 0
    @assert ecutrho >= 0
    @assert ecfixed >= 0
    @assert qcutz >= 0
    @assert q2sigma >= 0
    @assert length(starting_charge) <= ntyp
    @assert length(starting_magnetization) <= ntyp
    @assert length(Hubbard_U) <= ntyp
    @assert length(Hubbard_J0) <= ntyp
    @assert length(Hubbard_alpha) <= ntyp
    @assert length(Hubbard_beta) <= ntyp
    # @assert all(length(x) <= ntyp for x in Hubbard_J)
    @assert length(angle1) <= ntyp
    @assert length(angle2) <= ntyp
    @assert length(fixed_magnetization) <= 3
    @assert length(london_c6) <= ntyp
    @assert length(london_rvdw) <= ntyp
    @assert(
        exxdiv_treatment ∈
        ("gygi-baldereschi", "gygi-bald", "g-b", "vcut_ws", "vcut_spherical", "none")
    )
    @assert(
        !(x_gamma_extrapolation && exxdiv_treatment ∈ ("vcut_ws", "vcut_spherical")),
        "`x_gamma_extrapolation` cannot be used with `vcut`!"
    )
end # struct SystemNamelist

"""
    ElectronsNamelist <: Namelist

Represent the `ELECTRONS` namelist of `pw.x`.
"""
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
    @assert mixing_mode ∈ ("plain", "TF", "local-TF")
    @assert diagonalization ∈ ("david", "cg", "cg-serial", "david-serial", "ppcg")  # Different from docs
    @assert efield_phase ∈ ("read", "write", "none")
    @assert startingpot ∈ ("atomic", "file")
    @assert startingwfc ∈ ("atomic", "atomic+random", "random", "file")
end # struct ElectronsNamelist

"""
    IonsNamelist <: Namelist

Represent the `IONS` namelist of `pw.x`.

Input this namelist only if `calculation` is `"relax"`, `"md"`, `"vc-relax"`, or `"vc-md"`.
"""
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
        ion_dynamics ∈
        ("none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman")
    )
    @assert ion_positions ∈ ("default", "from_input")
    @assert pot_extrapolation ∈ ("none", "atomic", "first_order", "second_order")
    @assert wfc_extrapolation ∈ ("none", "first_order", "second_order")
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
        )
    )
    @assert tempw > 0
end # struct IonsNamelist

"""
    CellNamelist <: Namelist

Represent the `CELL` namelist of `pw.x`.

Input this namelist only if `calculation` is `"vc-relax"` or `"vc-md"`.
"""
@with_kw struct CellNamelist <: Namelist
    cell_dynamics::String = "none"
    press::Float64 = 0.0
    wmass::Float64 = 0.0
    cell_factor::Float64 = 0.0
    press_conv_thr::Float64 = 0.5
    cell_dofree::String = "all"
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1596-L1625.
    @assert cell_dynamics ∈ ("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w")
    @assert wmass >= 0
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
        )
    )
end # struct CellNamelist

# The following default values are picked from `<QE source>/PP/src/dos.f90`
"""
    DosNamelist <: Namelist

Represent the `DOS` namelist of `dos.x`.
"""
@with_kw struct DosNamelist <: Namelist
    prefix::String = "pwscf"
    outdir::String = "./"
    ngauss::Int = 0
    degauss::Float64 = 0.0
    Emin::Float64 = -1000000.0
    Emax::Float64 = 1000000.0
    DeltaE::Float64 = 0.01
    fildos::String = "$(prefix).dos"
    @assert ngauss ∈ (0, 1, -1, -99)
end # struct DosNamelist

# The following default values are picked from `<QE source>/PP/src/bands.f90`
"""
    BandsNamelist <: Namelist

Represent the `BANDS` namelist of `bands.x`.
"""
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
    @assert spin_component ∈ 1:2
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

"""
    Crystallography.BravaisLattice(nml::SystemNamelist)

Return a `BravaisLattice` from a `SystemNamelist`.
"""
function Crystallography.BravaisLattice(nml::SystemNamelist)
    ibrav = nml.ibrav
    return iszero(ibrav) ? nothing : BravaisLattice{nml.ibrav}()
end # function Crystallography.BravaisLattice

"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
function Crystallography.cellvolume(nml::SystemNamelist)
    iszero(nml.ibrav) && error("`ibrav` must be non-zero to calculate the cell volume!")
    return abs(det(BravaisLattice(nml)()))
end # function Crystallography.cellvolume

function Setters.make(::LensMaker{<:VerbositySetter,ControlNamelist})
    return @batchlens begin
        _.verbosity
        _.wf_collect
        _.tstress
        _.tprnfor
        _.disk_io
    end
end # function Setters.make
function Setters.make(::LensMaker{<:FiniteTemperatureSetter,SystemNamelist})
    return @batchlens begin
        _.occupations
        _.degauss
        _.smearing
    end
end # function Setters.make
function Setters.make(::LensMaker{<:CalculationSetter,ControlNamelist})
    return @lens _.calculation
end # function Setters.make

Setters.preset_values(::VerbositySetter{:high}, ::ControlNamelist) =
    ("high", true, true, true, "high")
Setters.preset_values(::VerbositySetter{:low}, ::ControlNamelist) =
    ("low", false, false, false, "low")
function Setters.preset_values(::FiniteTemperatureSetter{N}, ::SystemNamelist) where {N}
    return (
        "smearing",
        if N isa Real
            N
        else  # A quantity
            u = unit(N) |> upreferred
            if u == u"Ry"
                N
            elseif u == u"kg*m^2*s^-2"  # u"hartree", u"J", u"eV", etc..
                uconvert(u"Ry", N)
            elseif u == u"s^-1"  # u"Hz", u"THz", ...
                uconvert(u"Hz", N) / 6579683920502000.0 * 2
            elseif u == u"K"  # u"K", u"mK", u"μK", ...
                uconvert(u"K", N) / 315775.02480407 * 2
            elseif u == u"m^-1"  # u"m^-1", u"cm^-1", ...
                uconvert(u"m^-1", N) / 21947463.13632 * 2
            elseif u == u"kg"  # u"kg", u"g", ...
                uconvert(u"kg", N) / 4.8508702095432e-35 * 2
            else
                error("unknown unit given!")
            end |> ustrip
        end,
        "fermi-dirac",
    )
end # function Setters.preset_values
Setters.preset_values(::CalculationSetter{T}, ::ControlNamelist) where {T} =
    replace(string(T), "_" => "-")

end
