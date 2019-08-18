"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Compat: isnothing
using IterTools: fieldvalues
using Parameters: @with_kw

using QuantumESPRESSOBase: bravais_lattice
using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs

export PWscfInput, autogenerate_cell_parameters, namelists, cards

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
end  # struct PWscfInput

function autogenerate_cell_parameters(obj::PWscfInput)
    return reconstruct(
        obj,
        Dict(
            :system => reconstruct(obj.system, ibrav = 0),
            :cell_parameters => reconstruct(obj.cell_parameters, data = bravais_lattice(system))
        )
    )
end # function autogenerate_cell_parameters

filter_field_by_supertype(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWscfInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWscfInput) = filter_field_by_supertype(input, Card)

end
