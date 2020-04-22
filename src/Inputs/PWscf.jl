"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: isnothing, eachrow
using Crystallography: Bravais, Lattice, CellParameters, Cell
using Crystallography.Arithmetics: cellvolume
using Formatting: sprintf1
using Kaleido: @batchlens
using LinearAlgebra: det
using Parameters: @with_kw
using Pseudopotentials: pseudopot_format
using Setfield: PropertyLens, get, set, @lens, @set
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful
using UnitfulAtomic

using ..Inputs:
    InputEntry, Namelist, Input, entryname, Card, getoption, allowed_options, qestring
using ...Setters:
    AlatPressSetter,
    LensMaker,
    VerbositySetter,
    FiniteTemperatureSetter,
    CellParametersSetter,
    CalculationSetter,
    LensMaker

import Crystallography
import Crystallography.Arithmetics
import Pseudopotentials
import ..Inputs
import ...Setters

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    PWInput,
    optconvert

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
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
    nstep::UInt = 50
    iprint::UInt = 100000
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
    nberrycyc::UInt = 1
    lorbm::Bool = false
    lberry::Bool = false
    gdir::UInt = 1  # The QE default value is `0`!
    nppstr::UInt = 0
    lfcpopt::Bool = false
    gate::Bool = false
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1282-L1369.
    @assert calculation ∈ ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
    @assert verbosity ∈ ("high", "low", "debug", "medium", "default", "minimal")
    @assert restart_mode ∈ ("from_scratch", "restart")
    @assert iprint >= 1
    @assert disk_io ∈ ("high", "medium", "low", "none", "default")
    @assert dt >= 0
    @assert max_seconds >= 0
    @assert etot_conv_thr >= 0
    @assert forc_conv_thr >= 0
    @assert gdir ∈ 1:3
    @assert(
        !all((gate, tefield, !dipfield)),
        "`gate` cannot be used with `tefield` if dipole correction is not active!"
    )
    @assert(
        !all((gate, dipfield, !tefield)),
        "dipole correction is not active if `tefield = false`!"
    )
end # struct ControlNamelist

"""
    SystemNamelist <: Namelist

Represent the `SYSTEM` namelist of `pw.x`.
"""
@with_kw struct SystemNamelist <: Namelist
    ibrav::Int = 0  # The default value in QE's source code is -1
    celldm::Vector{Maybe{Float64}} = zeros(6)  # Must specify
    A::Float64 = 0.0
    B::Float64 = 0.0
    C::Float64 = 0.0
    cosAB::Float64 = 0.0
    cosAC::Float64 = 0.0
    cosBC::Float64 = 0.0
    nat::UInt = 0
    ntyp::UInt = 0
    nbnd::UInt = 0
    tot_charge::Float64 = 0.0
    starting_charge::Vector{Maybe{Float64}} = []
    tot_magnetization::Float64 = -1.0
    starting_magnetization::Vector{Maybe{Float64}} = []
    ecutwfc::Float64 = 0.0
    ecutrho::Float64 = 0.0
    ecutfock::Float64 = 0.0
    nr1::UInt = 0
    nr2::UInt = 0
    nr3::UInt = 0
    nr1s::UInt = 0
    nr2s::UInt = 0
    nr3s::UInt = 0
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
    nspin::UInt = 1
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
    nqx1::UInt = 1
    nqx2::UInt = 1
    nqx3::UInt = 1
    localization_thr::Float64 = 0.0  # This is only for QE 6.4
    lda_plus_u::Bool = false
    lda_plus_u_kind::UInt = 0
    Hubbard_U::Vector{Maybe{Float64}} = []
    Hubbard_J0::Vector{Maybe{Float64}} = []
    Hubbard_alpha::Vector{Maybe{Float64}} = []
    Hubbard_beta::Vector{Maybe{Float64}} = []
    # Hubbard_J::Vector{Vector{Maybe{Float64}}} = [zeros(ntyp)]  # The default value in QE's source code is just one 0.0
    starting_ns_eigenvalue::Float64 = -1.0  # It's actually a multidimensional array.
    U_projection_type::String = "atomic"
    edir::UInt = 1
    emaxpos::Float64 = 0.5
    eopreg::Float64 = 0.1
    eamp::Float64 = 0.001  # The default value in QE's source code is 0.0
    angle1::Vector{Maybe{Float64}} = []
    angle2::Vector{Maybe{Float64}} = []
    constrained_magnetization::String = "none"
    fixed_magnetization::Vector{Maybe{Float64}} = zeros(3)  # The default value in QE's source code is just one 0.0
    lambda::Float64 = 1.0
    report::UInt = 100
    lspinorb::Bool = false
    assume_isolated::String = "none"
    esm_bc::String = "pbc"
    esm_w::Float64 = 0.0
    esm_efield::Float64 = 0.0
    esm_nfit::UInt = 4
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
    space_group::UInt = 0
    uniqueb::Bool = false
    origin_choice::UInt = 1
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
        if ibrav == 0
            true  # Skip the check, cannot use `nothing`
        else
            celldm[1] != 0 || A != 0  # Cannot use `iszero` to compare now!
        end,
        "invalid lattice parameters (`celldm` $celldm or `A` $A)!"
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
    @assert(ntyp <= 10, "`ntyp` $ntyp is larger than 10!")
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
    @assert lda_plus_u_kind ∈ 0:1
    @assert edir ∈ 1:3
    @assert origin_choice ∈ 1:2
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
    electron_maxstep::UInt = 100
    scf_must_converge::Bool = true
    conv_thr::Float64 = 1e-6
    adaptive_thr::Bool = false
    conv_thr_init::Float64 = 1e-3
    conv_thr_multi::Float64 = 0.1
    mixing_mode::String = "plain"
    mixing_beta::Float64 = 0.7
    mixing_ndim::UInt = 8
    mixing_fixed_ns::UInt = 0
    diagonalization::String = "david"
    ortho_para::UInt = 0
    diago_thr_init::Float64 = 0.0
    diago_cg_maxiter::UInt = 20
    diago_david_ndim::UInt = 4
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
    nraise::UInt = 1
    refold_pos::Bool = false
    upscale::Float64 = 100.0
    bfgs_ndim::UInt = 1
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
    firstk::UInt = 0
    lastk::UInt = 10000000
    @assert spin_component ∈ 1:2
end # struct BandsNamelist

"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)

Represent each line of the `ATOMIC_SPECIES` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")

julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```
"""
struct AtomicSpecies
    "Label of the atom. Max total length cannot exceed 3 characters."
    atom::String
    """
    Mass of the atomic species in atomic unit.

    Used only when performing molecular dynamics (MD) run
    or structural optimization runs using damped MD.
    Not actually used in all other cases (but stored
    in data files, so phonon calculations will use
    these values unless other values are provided).
    """
    mass::Float64
    """
    File containing pseudopotential for this species.

    See also: [`pseudopot_format`](@ref)
    """
    pseudopot::String
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, pseudopot)
        @assert(length(atom) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        return new(string(atom), mass, pseudopot)
    end
end

"""
    pseudopot_format(data::AtomicSpecies)

Return the pseudopotential format of the `AtomicSpecies`.
"""
Pseudopotentials.pseudopot_format(data::AtomicSpecies) = pseudopot_format(data.pseudopot)

"""
    AtomicSpeciesCard <: Card

Represent the `ATOMIC_SPECIES` card in QE. It does not have an "option".
"""
struct AtomicSpeciesCard <: Card
    data::Vector{AtomicSpecies}
end
AtomicSpeciesCard(cell::Cell) = AtomicSpeciesCard(map(AtomicSpecies ∘ string, cell.numbers))

"""
    AtomicPosition(atom::Union{AbstractChar,String}, pos::Vector{Float64}[, if_pos::Vector{Int}])
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Represent each line of the `ATOMIC_POSITIONS` card in QE.

The `atom` field accepts at most 3 characters.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf

julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])

julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```
"""
struct AtomicPosition
    "Label of the atom as specified in `AtomicSpecies`."
    atom::String
    "Atomic positions. A three-element vector of floats."
    pos::SVector{3,Float64}
    """
    Component `i` of the force for this atom is multiplied by `if_pos(i)`,
    which must be either `0` or `1`.  Used to keep selected atoms and/or
    selected components fixed in MD dynamics or structural optimization run.

    With `crystal_sg` atomic coordinates the constraints are copied in all equivalent
    atoms.
    """
    if_pos::SVector{3,Bool}
    function AtomicPosition(atom::Union{AbstractChar,AbstractString}, pos, if_pos)
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom), pos, if_pos)
    end
end
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicPosition, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)

"""
    AtomicPositionsCard <: Card

Represent the `ATOMIC_POSITIONS` card in QE.

# Arguments
- `data::AbstractVector{AtomicPosition}`: A vector containing `AtomicPosition`s.
- `option::String="alat"`: allowed values are: "alat", "bohr", "angstrom", "crystal", and "crystal_sg".
"""
struct AtomicPositionsCard <: Card
    data::Vector{AtomicPosition}
    option::String
    function AtomicPositionsCard(data, option = "alat")
        @assert option ∈ allowed_options(AtomicPositionsCard)
        return new(data, option)
    end
end
AtomicPositionsCard(cell::Cell, option) = AtomicPositionsCard(
    [
        AtomicPosition(string(atom), pos)
        for (atom, pos) in zip(cell.numbers, cell.positions)
    ],
    option,
)
# Introudce mutual constructors since they share the same atoms.

function validate(x::AtomicSpeciesCard, y::AtomicPositionsCard)
    lens = @lens _.data.atom
    @assert(
        isempty(symdiff(map(Base.Fix2(get, lens) ∘ unique, (x, y)))),
        "labels of the atoms are different in `ATOMIC_SPECIES` and `ATOMIC_POSITIONS` card!",
    )
end # function validate
validate(y::AtomicPositionsCard, x::AtomicSpeciesCard) = validate(x, y)

"Represent the abstraction of `CELL_PARAMETERS` and `REF_CELL_PARAMETERS` cards in QE."
abstract type AbstractCellParametersCard <: Card end

"""
    CellParametersCard{T<:Real} <: AbstractCellParametersCard
    CellParametersCard(data::AbstractMatrix, option::String)

Represent the `CELL_PARAMETERS` cards in `PWscf` and `CP` packages.
"""
struct CellParametersCard{T<:Real} <: AbstractCellParametersCard
    data::SMatrix{3,3,T}
    option::String
    function CellParametersCard{T}(data, option = "alat") where {T<:Real}
        @assert option ∈ allowed_options(CellParametersCard)
        return new(data, option)
    end
end
CellParametersCard(data::AbstractMatrix{T}, option = "alat") where {T} =
    CellParametersCard{T}(data, option)
CellParametersCard(lattice::Lattice{T}, option) where {T} =
    CellParametersCard(convert(Matrix{T}, lattice), option)
CellParametersCard(cell::Cell, option) = CellParametersCard(cell.lattice, option)

struct AtomicForce
    atom::String
    force::SVector{3,Float64}
    function AtomicForce(atom::Union{AbstractChar,AbstractString}, force)
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom), force)
    end
end

struct AtomicForcesCard <: Card
    data::Vector{AtomicForce}
end

"""
    optconvert(new_option::AbstractString, card::AbstractCellParametersCard)

Convert the option of an `AbstractCellParametersCard` from "bohr" to "angstrom", or its reverse.

!!! warning
    It does not support conversion between `"alat"` and the others.
"""
function optconvert(new_option::AbstractString, card::AbstractCellParametersCard)
    old_option = getoption(card)
    if new_option == old_option
        return card  # No conversion is needed
    else
        return typeof(card)(
            if (old_option => new_option) == ("bohr" => "angstrom")
                @. ustrip(u"angstrom", card.data * u"bohr")
            elseif (old_option => new_option) == ("angstrom" => "bohr")
                @. ustrip(u"bohr", card.data * u"angstrom")
            else
                error("unknown conversion rule $(old_option => new_option)!")
            end,
        )
    end
end # function optconvert

"""
    MonkhorstPackGrid(grid, offsets)

Represent the Monkhorst--Pack grid.

# Arguments
- `grid`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `offsets`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    grid::SVector{3,UInt}
    offsets::SVector{3,Bool}
end

"Represent the centre of the Brillouin zone (commonly marked as the Γ point)."
struct GammaPoint end

"""
    SpecialKPoint(coord, weight)

Represent a special point of the 3D Brillouin zone. Each of them has a weight.
"""
struct SpecialKPoint <: FieldVector{4,Float64}
    x::Float64
    y::Float64
    z::Float64
    w::Float64
end
SpecialKPoint(::GammaPoint) = SpecialKPoint(0.0, 0.0, 0.0, 1.0)

"""
    struct KPointsCard{<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}} <: Card

Represent the `K_POINTS` card in QE.

# Arguments
- `data::Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}`: A Γ point, a Monkhorst--Pack grid or a vector containing `SpecialKPoint`s.
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
"""
struct KPointsCard{A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}} <:
       Card
    data::A
    option::String
    function KPointsCard{A}(
        data,
        option,
    ) where {A<:Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}}
        @assert option ∈ allowed_options(KPointsCard)
        @assert if option == "automatic"
            typeof(data) <: MonkhorstPackGrid
        elseif option == "gamma"
            typeof(data) <: GammaPoint
        else  # option ∈ ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            eltype(data) <: SpecialKPoint
        end
        return new(data, option)
    end
end
KPointsCard(data::A, option) where {A} = KPointsCard{A}(data, option)
KPointsCard(data::AbstractVector{SpecialKPoint}) = KPointsCard(data, "tpiba")
KPointsCard(data::GammaPoint) = KPointsCard(data, "gamma")
KPointsCard(data::MonkhorstPackGrid) = KPointsCard(data, "automatic")

"""
    Crystallography.Bravais(nml::SystemNamelist)

Return a `Bravais` from a `SystemNamelist`.
"""
Crystallography.Bravais(nml::SystemNamelist) = Bravais(nml.ibrav)

"""
    Crystallography.Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
function Crystallography.Lattice(nml::SystemNamelist)
    b = Bravais(nml)
    return Lattice(b, CellParameters(nml.celldm...))
end # function Crystallography.Lattice

"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function Arithmetics.cellvolume(card::AbstractCellParametersCard)
    option = getoption(card)
    if option == "bohr"
        abs(det(card.data))
    elseif option == "angstrom"
        ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
        error("information not enough! Parameter `celldm[1]` needed!")
    end
end # function Arithmetics.cellvolume
"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
Arithmetics.cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))

function Setters.make(::LensMaker{CellParametersSetter})
    return @batchlens begin
        _.cell_parameters
        _.system.ibrav
        _.system.celldm
    end
end # function Setters.make

function Setters.preset_values(::CellParametersSetter, template)
    # !isnothing(template.cell_parameters) && return template
    system = template.system
    return (
        CellParametersCard(
            Lattice(Bravais(system), CellParameters(template.celldm...)),
            "alat",
        ),
        0,
        [system.celldm[1]],
    )
end # function Setters.preset_values

Inputs.getoption(::AtomicSpeciesCard) = nothing

Inputs.allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
Inputs.allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
Inputs.allowed_options(::Type{<:KPointsCard}) = (
    "tpiba",
    "automatic",
    "crystal",
    "gamma",
    "tpiba_b",
    "crystal_b",
    "tpiba_c",
    "crystal_c",
)

Inputs.titleof(::Type{ControlNamelist}) = "CONTROL"
Inputs.titleof(::Type{SystemNamelist}) = "SYSTEM"
Inputs.titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
Inputs.titleof(::Type{IonsNamelist}) = "IONS"
Inputs.titleof(::Type{CellNamelist}) = "CELL"
Inputs.titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
Inputs.titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
Inputs.titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"
Inputs.titleof(::Type{<:KPointsCard}) = "K_POINTS"

function Inputs.qestring(data::AtomicSpecies; delim = ' ', numfmt = "%14.9f", args...)
    return join(
        (sprintf1("%3s", data.atom), sprintf1(numfmt, data.mass), data.pseudopot),
        delim,
    )
end
function Inputs.qestring(
    card::AtomicSpeciesCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%20.10f",
    newline = '\n',
)
    # Using generator expressions in `join` is faster than using `Vector`s.
    return "ATOMIC_SPECIES" *
           newline *
           join(
               (
                   indent * qestring(x; delim = delim, numfmt = numfmt)
                   for x in unique(card.data)
               ),
               newline,
           )
end
function Inputs.qestring(data::AtomicPosition; delim = ' ', numfmt = "%14.9f", args...)
    f(x) = x ? "" : "0"
    return join(
        [
            sprintf1("%3s", data.atom)
            map(x -> sprintf1(numfmt, x), data.pos)
            map(f, data.if_pos)
        ],
        delim,
    )
end
function Inputs.qestring(
    card::AtomicPositionsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    return "ATOMIC_POSITIONS { $(getoption(card)) }" *
           newline *
           join(
               (indent * qestring(x; delim = delim, numfmt = numfmt) for x in card.data),
               newline,
           )
end
function Inputs.qestring(
    card::CellParametersCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    it = (
        indent * join((sprintf1(numfmt, x) for x in row), delim) for
        row in eachrow(card.data)
    )
    return "CELL_PARAMETERS { $(getoption(card)) }" * newline * join(it, newline)
end
Inputs.qestring(data::GammaPoint) = ""
Inputs.qestring(data::MonkhorstPackGrid; delim = ' ', numfmt = "%5d", args...) =
    join(map(x -> sprintf1(numfmt, x), [data.grid; data.offsets]), delim)
Inputs.qestring(data::SpecialKPoint; delim = ' ', numfmt = "%14.9f", args...) =
    join(map(x -> sprintf1(numfmt, x), collect(data)), delim)
function Inputs.qestring(
    card::KPointsCard;
    indent = ' '^4,
    delim = ' ',
    numfmt = "%14.9f",
    newline = '\n',
)
    content = "K_POINTS { $(card.option) }" * newline
    if getoption(card) in ("gamma", "automatic")
        content *= indent * qestring(card.data)
    else  # ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= string(length(card.data), newline)
        content *= join(
            (indent * qestring(x; delim = delim, numfmt = numfmt) for x in card.data),
            newline,
        )
    end
    return content
end

function Base.setproperty!(value::AtomicSpecies, name::Symbol, x)
    if name == :atom
        @assert(length(x) <= 3, "Max total length of `atom` cannot exceed 3 characters!")
        x = string(x)  # An `if` statement is more expensive than directly setting a string
    end
    setfield!(value, name, x)  # FIXME: It is now immutable!
end # function Base.setproperty!
function Base.setproperty!(value::AtomicPosition, name::Symbol, x)
    x = if name == :atom
        @assert(length(x) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        x = string(x)  # Make sure it is a `String`
    elseif name ∈ (:pos, :if_pos) && x isa AbstractVector
        SVector{3}(x)
    end
    setfield!(value, name, x)  # FIXME: It is now immutable!
end # function Base.setproperty!

"""
    PWInput <: Input
    PWInput(control, system, electrons, ions, cell, atomic_species, atomic_positions, k_points, cell_parameters)

Construct a `PWInput` which represents the input of program `pw.x`.

# Arguments
- `control::ControlNamelist=ControlNamelist()`: the `CONTROL` namelist of the input. Optional.
- `system::SystemNamelist=SystemNamelist()`: the `SYSTEM` namelist of the input. Optional.
- `electrons::ElectronsNamelist=ElectronsNamelist()`: the `ELECTRONS` namelist of the input. Optional.
- `ions::IonsNamelist=IonsNamelist()`: the `IONS` namelist of the input. Optional.
- `cell::CellNamelist=CellNamelist()`: the `CELL` namelist of the input. Optional.
- `atomic_species::AtomicSpeciesCard`: the `ATOMIC_SPECIES` card of the input. Must be provided explicitly.
- `atomic_positions::AtomicPositionsCard`: the `ATOMIC_POSITIONS` card of the input. Must be provided explicitly.
- `k_points::KPointsCard`: the `K_POINTS` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
@with_kw struct PWInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard} = nothing
    constraints::Union{Union{Nothing,Float64}} = nothing
    occupations::Union{Nothing,Float64} = nothing
    atomic_forces::Union{Nothing,AtomicForcesCard} = nothing
    @assert(
        !(isnothing(cell_parameters) && system.ibrav == 0),
        "Cannot specify an empty `cell_parameters` with `ibrav = 0`!"
    )
end # struct PWInput
PWInput(args...) =
    PWInput(; Dict(zip(map(arg -> entryname(typeof(arg), PWInput), args), args))...)  # See https://discourse.julialang.org/t/construct-an-immutable-type-from-a-dict/26709/6

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

function Setters.make(::LensMaker{AlatPressSetter,PWInput})
    return @batchlens begin
        _.system.celldm  # Get the `template`'s `system.celldm` value
        _.cell.press     # Get the `template`'s `cell.press` value
        _.cell_parameters.option
    end
end # function Setters.make
# function Setters.upgrade(lm::LensMaker{S,T}, ::Type{PWInput}) where {S,T<:InputEntry}
# return PropertyLens{entryname(T)}() ∘ make(lm)
# end # function Setters.upgrade

end
