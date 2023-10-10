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
    fcp::FcpNamelist
    rism::RismNamelist
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Maybe{CellParametersCard}
    occupations::Maybe{OccupationsCard}
    constraints::Maybe{ConstraintsCard}
    atomic_velocities::Maybe{AtomicVelocitiesCard}
    atomic_forces::Maybe{AtomicForcesCard}
    additional_k_points::Maybe{AdditionalKPointsCard}
    solvents::Maybe{SolventsCard}
    hubbard::Maybe{HubbardCard}
end
function PWInput(;
    control=ControlNamelist(),
    system,
    electrons=ElectronsNamelist(),
    ions=IonsNamelist(),
    cell=CellNamelist(),
    fcp=FcpNamelist(),
    rism=RismNamelist(),
    atomic_species,
    atomic_positions,
    k_points,
    cell_parameters=nothing,
    occupations=nothing,
    constraints=nothing,
    atomic_velocities=nothing,
    atomic_forces=nothing,
    additional_k_points=nothing,
    solvents=nothing,
    hubbard=nothing,
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
        fcp,
        rism,
        atomic_species,
        atomic_positions,
        k_points,
        cell_parameters,
        occupations,
        constraints,
        atomic_velocities,
        atomic_forces,
        additional_k_points,
        solvents,
        hubbard,
    )
end

exitfile(input::PWInput) = exitfile(input.control)

mkexitfile(input::PWInput) = mkexitfile(input.control)

"""
    isrequired(nml::Namelist)
    isrequired(card::Card)

Test whether a `Namelist` or a `Card` is required in `PWInput`.
"""
isrequired(nml::Namelist) = nml isa Union{ControlNamelist,SystemNamelist,ElectronsNamelist}
isrequired(card::Card) = card isa Union{AtomicSpeciesCard,AtomicPositionsCard,KPointsCard}

"""
    isoptional(nml::Namelist)
    isoptional(card::Card)

Test whether a `Namelist` or a `Card` is optional in `PWInput`.
"""
isoptional(nml::Namelist) =
    nml isa Union{IonsNamelist,CellNamelist,FcpNamelist,RismNamelist}
isoptional(card::Card) =
    card isa Union{
        CellParametersCard,
        OccupationsCard,
        ConstraintsCard,
        AtomicVelocitiesCard,
        AtomicForcesCard,
        AdditionalKPointsCard,
        SolventsCard,
        HubbardCard,
    }

"""
    eachpotential(input::PWInput)

Iterate the pseudopotentials in a `PWInput`.
"""
eachpotential(input::PWInput) = eachpotential(input.atomic_species)

"""
    getpseudodir(input::PWInput)

Get the directory storing the pseudopotential files.
"""
getpseudodir(input::PWInput) = getpseudodir(input.control)

getxmldir(input::PWInput) = getxmldir(input.control)
