"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using IterTools: fieldvalues
using Parameters: @with_kw

using QuantumESPRESSOBase: bravais_lattice
using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs

export PWscfInput, namelists, cards

@with_kw struct PWscfInput <: AbstractInput
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist = SystemNamelist()
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::Union{Nothing,CellParametersCard} = nothing
    function PWscfInput(
        control::ControlNamelist,
        system::SystemNamelist,
        electrons::ElectronsNamelist,
        ions::IonsNamelist,
        cell::CellNamelist,
        atomic_species::AtomicSpeciesCard,
        atomic_positions::AtomicPositionsCard,
        k_points::KPointsCard,
        cell_parameters::Nothing
    )
        system.ibrav == 0 && error("Cannot specify `ibrav = 0` with an empty `cell_parameters`!")
        new(control, system, electrons, ions, cell, atomic_species, atomic_positions, k_points, bravais_lattice(system))
    end
end  # struct PWscfInput

filter_field_by_supertype(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWscfInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWscfInput) = filter_field_by_supertype(input, Card)

end
