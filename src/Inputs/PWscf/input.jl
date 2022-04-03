import ..Inputs:
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards

export PWInput,
    exitfile,
    mkexitfile,
    allnamelists,
    allcards,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards

"""
    PWInput <: QuantumESPRESSOInput
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
- `k_points::AbstractKPointsCard`: the `K_POINTS` card of the input. Must be provided explicitly.
- `cell_parameters::Union{Nothing,CellParametersCard}`: the `CELL_PARAMETERS` card of the input. Must be either `nothing` or a `CellParametersCard`.
"""
struct PWInput <: QuantumESPRESSOInput
    control::ControlNamelist
    system::SystemNamelist
    electrons::ElectronsNamelist
    ions::IonsNamelist
    cell::CellNamelist
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard}
    constraints::Union{Union{Nothing,Float64}}
    occupations::Union{Nothing,Float64}
    atomic_forces::Union{Nothing,AtomicForcesCard}
end # struct PWInput
function PWInput(;
    control = ControlNamelist(),
    system,
    electrons = ElectronsNamelist(),
    ions = IonsNamelist(),
    cell = CellNamelist(),
    atomic_species,
    atomic_positions,
    k_points,
    cell_parameters = nothing,
    constraints = nothing,
    occupations = nothing,
    atomic_forces = nothing,
)
    @assert !isnothing(cell_parameters) || system.ibrav != 0 "`cell_parameters` is empty with `ibrav = 0`!"
    foreach(atomic_species.data) do datum
        path = joinpath(control.pseudo_dir, datum.pseudopot)
        if !isfile(path)
            @warn "pseudopotential file \"$path\" does not exist!"
        end
    end
    return PWInput(
        control,
        system,
        electrons,
        ions,
        cell,
        atomic_species,
        atomic_positions,
        k_points,
        cell_parameters,
        constraints,
        occupations,
        atomic_forces,
    )
end

@batteries PWInput eq = true hash = true

exitfile(template::PWInput) = abspath(
    expanduser(joinpath(template.control.outdir, template.control.prefix * ".EXIT")),
)
function mkexitfile(template::PWInput)
    path = exitfile(template)
    mkpath(dirname(path))
    return touch(path)
end

"""
    allnamelists(input::PWInput)

Return an iterator of all `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
allnamelists(input::PWInput) =
    (getfield(input, f) for f in (:control, :system, :electrons, :ions, :cell))

"""
    allcards(input::PWInput)

Get all `Card`s from a `PWInput`.
"""
allcards(input::PWInput) = (
    getfield(input, f) for f in (
        :atomic_species,
        :atomic_positions,
        :k_points,
        :cell_parameters,
        :constraints,
        :occupations,
        :atomic_forces,
    )
)

"""
    required_namelists(input::PWInput)

Return an iterator of required `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
required_namelists(input::PWInput) = (getfield(input, f) for f in (:control, :system, :electrons))

"""
    optional_namelists(input::PWInput)

Return an iterator of optional `Namelist`s from a `PWInput`. You may want to `collect` them.
"""
optional_namelists(input::PWInput) = (getfield(input, f) for f in (:ions, :cell))

"""
    required_cards(input::PWInput)

Return an iterator of required `Card`s from a `PWInput`. You may want to `collect` them.
"""
required_cards(input::PWInput) =
    (getfield(input, f) for f in (:atomic_species, :atomic_positions, :k_points))

"""
    optional_cards(input::PWInput)

Return an iterator of optional `Card`s from a `PWInput`. You may want to `collect` them.
"""
optional_cards(input::PWInput) =
    (getfield(input, f) for f in (:cell_parameters, :constraints, :occupations, :atomic_forces))

"""
    getpotentials(input::PWInput)

Get the pseudopotential names from a `PWInput`.
"""
getpotentials(input::PWInput) = getpotentials(input.atomic_species)

"""
    getpseudodir(input::PWInput)

Get the directory storing the pseudopotential files.
"""
getpseudodir(input::PWInput) = getpseudodir(input.control)
