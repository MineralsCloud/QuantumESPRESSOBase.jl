module CP

using Compat: isnothing
using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists.CP:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    PressAiNamelist,
    WannierNamelist
using QuantumESPRESSOBase.Cards.CP:
    AtomicSpeciesCard,
    AtomicPositionsCard,
    AtomicVelocitiesCard,
    CellParametersCard,
    RefCellParametersCard,
    AtomicForcesCard
using QuantumESPRESSOBase.Inputs: QuantumESPRESSOInput

export CPInput

"""
    CPInput <: QuantumESPRESSOInput
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
@with_kw struct CPInput <: QuantumESPRESSOInput
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
