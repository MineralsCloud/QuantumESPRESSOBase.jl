export PWInput, isrequired, isoptional

"""
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
@struct_hash_equal struct PWInput <: QuantumESPRESSOInput
    control::ControlNamelist
    system::SystemNamelist
    electrons::ElectronsNamelist
    ions::IonsNamelist
    cell::CellNamelist
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Maybe{CellParametersCard}
    constraints::Maybe{Float64}
    occupations::Maybe{Float64}
    atomic_forces::Maybe{AtomicForcesCard}
end
function PWInput(;
    control=ControlNamelist(),
    system,
    electrons=ElectronsNamelist(),
    ions=IonsNamelist(),
    cell=CellNamelist(),
    atomic_species,
    atomic_positions,
    k_points,
    cell_parameters=nothing,
    constraints=nothing,
    occupations=nothing,
    atomic_forces=nothing,
)
    @assert !isnothing(cell_parameters) || system.ibrav != 0 "`cell_parameters` is empty with `ibrav = 0`!"
    foreach(eachatom(atomic_species)) do atom
        path = joinpath(expanduser(control.pseudo_dir), atom.pseudopot)
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

exitfile(input::PWInput) = exitfile(input.control)

mkexitfile(input::PWInput) = mkexitfile(input.control)

isrequired(nml::Namelist) = nml isa Union{ControlNamelist,SystemNamelist,ElectronsNamelist}
isrequired(card::Card) = card isa Union{AtomicSpeciesCard,AtomicPositionsCard,KPointsCard}

isoptional(nml::Namelist) = nml isa Union{IonsNamelist,CellNamelist}
isoptional(card::Card) = card isa Union{CellParametersCard,AtomicForcesCard}

"""
    getpotentials(input::PWInput)

Get the pseudopotential names from a `PWInput`.
"""
listpotentials(input::PWInput) = listpotentials(input.atomic_species)

"""
    getpseudodir(input::PWInput)

Get the directory storing the pseudopotential files.
"""
getpseudodir(input::PWInput) = getpseudodir(input.control)

getxmldir(input::PWInput) = getxmldir(input.control)
