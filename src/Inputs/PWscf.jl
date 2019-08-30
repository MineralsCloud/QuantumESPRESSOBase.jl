"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: isnothing
using Parameters: @with_kw, reconstruct

using QuantumESPRESSOBase: bravais_lattice
using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs

export PWscfInput, autofill_cell_parameters, namelists, cards, compulsory_namelists, compulsory_cards

"""
    PWscfInput(control, system, electrons, ions, cell, atomic_species, atomic_positions, k_points, cell_parameters)

Construct a `PWscfInput` which represents the input of program `pw.x`.

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
@with_kw struct PWscfInput <: AbstractInput
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard}
    @assert !(isnothing(cell_parameters) && system.ibrav == 0) "Cannot specify `ibrav = 0` with an empty `cell_parameters`!"
end # struct PWscfInput

"""
    autofill_cell_parameters(obj::PWscfInput)

Generate automatically a `CellParametersCard` for a `PWscfInput` if its `cell_parameters` field is `nothing`.

Sometimes the `ibrav` field of a `PWscfInput` is not `0`, with its `cell_parameters` field to be empty.
But there are cases we want to write its `CellParametersCard` explicitly. This function will take a `PWscfInput` described
above and generate a new `PWscfInput` with its `ibrav = 0` and `cell_parameters` not empty.
"""
function autofill_cell_parameters(obj::PWscfInput)
    return reconstruct(
        obj,
        Dict(
            :system => reconstruct(obj.system, ibrav = 0),
            :cell_parameters => reconstruct(obj.cell_parameters, data = bravais_lattice(system))
        )
    )
end # function autofill_cell_parameters

"""
    filter_field_by_supertype(obj, ::Type)

A helper function to implement `namelists` and `cards`. It should not be exported.
"""
filter_field_by_supertype(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

"""
    namelists(input::PWscfInput)

Return a vector of `Namelist`s of a `PWscfInput`.
"""
namelists(input::PWscfInput) = filter_field_by_supertype(input, Namelist)

"""
    cards(input::PWscfInput)

Return a vector of `Card`s of a `PWscfInput`.
"""
cards(input::PWscfInput) = filter_field_by_supertype(input, Card)

"""
    compulsory_namelists(input::PWscfInput)

Return a vector of compulsory `Namelist`s of a `PWscfInput`.

The compulsory `Namelist`s of a `PWscfInput` are `ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`.
"""
compulsory_namelists(input::PWscfInput) = [getfield(input, x) for x in (:control, :system, :electrons)]

"""
    compulsory_cards(input::PWscfInput)

Return a vector of compulsory `Card`s of a `PWscfInput`.

The compulsory `Card`s of a `PWscfInput` are `AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`.
"""
compulsory_cards(input::PWscfInput) = [getfield(input, x) for x in (:atomic_species, :atomic_positions, :k_points)]

end
