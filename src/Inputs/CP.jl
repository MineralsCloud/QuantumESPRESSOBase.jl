module CP

using Compat: isnothing, eachrow
using Crystallography: Bravais, Lattice, CellParameters, Cell
using Formatting: sprintf1
using Kaleido: @batchlens
using LinearAlgebra: det
using Parameters: @with_kw
using Pseudopotentials: pseudopot_format
using Setfield: get, set, @lens, @set
using StaticArrays: SVector, SMatrix, FieldVector
using Unitful
using UnitfulAtomic

using ..Inputs: Namelist, Input, Card, getoption, allowed_options, qestring
using ...Setters: VerbositySetter, CalculationSetter, LensMaker

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
    PressAiNamelist,
    WannierNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    RefCellParametersCard,
    AtomicVelocity,
    AtomicVelocitiesCard,
    AtomicForce,
    AtomicForcesCard,
    CPInput,
    optconvert

# From https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4
const Maybe{T} = Union{T,Nothing}

# The following default values are picked from `<QE source>/Modules/read_namelists.f90`
"""
    ControlNamelist <: Namelist

Represent the `CONTROL` namelist of `cp.x`.
"""
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
        calculation ∈
        ("cp", "scf", "nscf", "relax", "vc-relax", "vc-cp", "cp-wf", "vc-cp-wf")
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

"""
    SystemNamelist <: Namelist

Represent the `SYSTEM` namelist of `cp.x`.
"""
@with_kw struct SystemNamelist <: Namelist
    ibrav::Int = -1
    celldm::Vector{Maybe{Float64}} = zeros(6)  # Must specify
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
    @assert ibrav ∈ union(0:1:14, (-3, -5, -9, -12, -13))
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
        elseif ibrav ∈ (4, 6, 7, 8, 9, -9, 10, 11)
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
    @assert nspin ∈ (1, 2, 4)
    @assert ecutwfc >= 0
    @assert ecutrho >= 0
    @assert(iszero(degauss), "`degauss` is not used in CP!")
    @assert ecfixed >= 0
    @assert qcutz >= 0
    @assert q2sigma >= 0
end # struct SystemNamelist

"""
    ElectronsNamelist <: Namelist

Represent the `ELECTRONS` namelist of `cp.x`.
"""
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

"""
    IonsNamelist <: Namelist

Represent the `IONS` namelist of `cp.x`.

Input this namelist only if `calculation` is `"cp"`, `"relax"`, `"vc-relax"`, `"vc-cp"`, `"cp-wf"`, or `"vc-cp-wf"`.
"""
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

"""
    CellNamelist <: Namelist

Represent the `CELL` namelist of `cp.x`.

Input this namelist only if `calculation` is `"vc-relax"`, `"vc-cp"`, or `"vc-cp-wf"`.
"""
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

"""
    PressAiNamelist <: Namelist

Represent the `PRESS_AI` namelist of `cp.x`.

Input this namelist only when `tabps` is `true`.
"""
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

"""
    WannierNamelist <: Namelist

Represent the `WANNIER` namelist of `cp.x`.

Input this namelist only if `calculation` is `"cp-wf"` or `"vc-cp-wf"`.
"""
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

Inputs.getoption(::AtomicSpeciesCard) = nothing

Inputs.allowed_options(::Type{<:AtomicPositionsCard}) =
    ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
Inputs.allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")

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

Inputs.titleof(::Type{ControlNamelist}) = "CONTROL"
Inputs.titleof(::Type{SystemNamelist}) = "SYSTEM"
Inputs.titleof(::Type{ElectronsNamelist}) = "ELECTRONS"
Inputs.titleof(::Type{IonsNamelist}) = "IONS"
Inputs.titleof(::Type{CellNamelist}) = "CELL"
Inputs.titleof(::Type{AtomicSpeciesCard}) = "ATOMIC_SPECIES"
Inputs.titleof(::Type{AtomicPositionsCard}) = "ATOMIC_POSITIONS"
Inputs.titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"

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

function Setters.make(::LensMaker{VerbositySetter,ControlNamelist})
    return @batchlens begin
        _.verbosity
        _.wf_collect
        _.tstress
        _.tprnfor
        _.saverho
        _.disk_io
    end
end # function Setters.make
function Setters.make(::LensMaker{<:CalculationSetter,ControlNamelist})
    return @lens _.calculation
end # function Setters.make

Setters.preset_values(::VerbositySetter{:high}, ::ControlNamelist) =
    ("high", true, true, true, true, "high")
Setters.preset_values(::VerbositySetter{:low}, ::ControlNamelist) =
    ("low", false, false, false, false, "default")
Setters.preset_values(::CalculationSetter{T}, ::ControlNamelist) where {T} =
    replace(string(T), "_" => "-")

"""
    AtomicVelocity(atom::Union{AbstractChar,String}, velocity::Vector{Float64})
    AtomicVelocity(x::AtomicPosition, velocity)

Represent each line of the `ATOMIC_VELOCITIES` card in QE's `CP` package.

The `atom` field accepts at most 3 characters.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.CP

julia> AtomicVelocity("H", [0.140374e-04, -0.333683e-04, 0.231834e-04])
AtomicVelocity("H", [1.40374e-5, -3.33683e-5, 2.31834e-5])
```
"""
struct AtomicVelocity
    atom::String
    velocity::SVector{3,Float64}
    function AtomicVelocity(atom::Union{AbstractChar,AbstractString}, velocity)
        @assert(length(atom) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        return new(string(atom), velocity)
    end
end
AtomicVelocity(x::AtomicSpecies, velocity) = AtomicVelocity(x.atom, velocity)
AtomicVelocity(x::AtomicPosition, velocity) = AtomicVelocity(x.atom, velocity)
# Introudce mutual constructors since they share the same atoms.
"""
    AtomicSpecies(x::AtomicVelocity, mass, pseudopot)

Construct an `AtomicSpecies` from an `AtomicVelocity` instance.
"""
AtomicSpecies(x::AtomicVelocity, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)
"""
    AtomicPosition(x::AtomicVelocity, pos, if_pos)

Construct an `AtomicPosition` from an `AtomicVelocity` instance.
"""
AtomicPosition(x::AtomicVelocity, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)

"""
    AtomicVelocitiesCard <: Card

Represent the `ATOMIC_VELOCITIES` card in QE's `CP` package which does not have an "option".
"""
struct AtomicVelocitiesCard <: Card
    data::Vector{AtomicVelocity}
end

struct RefCellParametersCard{T<:Real} <: AbstractCellParametersCard
    data::SMatrix{3,3,T}
    option::String
    function RefCellParametersCard{T}(data, option = "bohr") where {T<:Real}
        @assert option ∈ allowed_options(RefCellParametersCard)
        return new(data, option)
    end
end
RefCellParametersCard(data::AbstractMatrix{T}, option = "bohr") where {T} =
    RefCellParametersCard{T}(data, option)

Inputs.getoption(::AtomicVelocitiesCard) = "a.u"
Inputs.getoption(::AtomicForcesCard) = nothing

Inputs.allowed_options(::Type{<:AtomicVelocity}) = ("a.u",)
Inputs.allowed_options(::Type{<:RefCellParametersCard}) = ("bohr", "angstrom")

function Base.setproperty!(value::AtomicVelocity, name::Symbol, x)
    x = if name == :atom
        @assert(length(x) <= 3, "the max length of `atom` cannot exceed 3 characters!")
        x = string(x)  # Make sure it is a `String`
    elseif name == :velocity && x isa AbstractVector
        SVector{3}(x)
    end
    setfield!(value, name, x)
end # function Base.setproperty!

"""
    CPInput <: Input
    CPInput(control, system, electrons, ions, cell, atomic_species, atomic_positions, k_points, cell_parameters)

Construct a `PWInput` which represents the input of program `pw.x`.

# Arguments
- `control::ControlNamelist=ControlNamelist()`: the `CONTROL` namelist of the input. Optional.
- `system::SystemNamelist=SystemNamelist()`: the `SYSTEM` namelist of the input. Optional.
- `electrons::ElectronsNamelist=ElectronsNamelist()`: the `ELECTRONS` namelist of the input. Optional.
- `ions::IonsNamelist=IonsNamelist()`: the `IONS` namelist of the input. Optional.
- `cell::CellNamelist=CellNamelist()`: the `CELL` namelist of the input. Optional.
- `press_ai::PressAiNamelist=PressAiNamelist()`: the `PRESS_AI` namelist of the input. Optional.
- `wannier::WannierNamelist=WannierNamelist()`: the `WANNIER` namelist of the input. Optional.
- `atomic_species::AtomicSpeciesCard`: the `ATOMIC_SPECIES` card of the input. Must be provided explicitly.
- `atomic_positions::AtomicPositionsCard`: the `ATOMIC_POSITIONS` card of the input. Must be provided explicitly.
- `atomic_velocities::AtomicVelocitiesCard`: the `ATOMIC_VELOCITIES` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
- `ref_cell_parameters::Union{Nothing,RefCellParametersCard}`: the `REF_CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
@with_kw struct CPInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    press_ai::PressAiNamelist = PressAiNamelist()
    wannier::WannierNamelist = WannierNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    atomic_velocities::AtomicVelocitiesCard
    cell_parameters::Union{Nothing,CellParametersCard} = nothing
    ref_cell_parameters::Union{Nothing,RefCellParametersCard} = nothing
    constraints::Union{Nothing,Float64} = nothing
    occupations::Union{Nothing,Float64} = nothing
    atomic_forces::Union{Nothing,AtomicForcesCard} = nothing
    plot_wannier::Union{Nothing,Float64} = nothing
    autopilot::Union{Nothing,Float64} = nothing
    @assert(
        !(isnothing(cell_parameters) && system.ibrav == 0),
        "Cannot specify `ibrav = 0` with an empty `cell_parameters`!"
    )
end # struct CPInput

end
