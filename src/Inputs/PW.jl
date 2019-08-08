"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using IterTools: fieldvalues
using Parameters: @with_kw

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PW
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PW
using QuantumESPRESSOBase.QuantumESPRESSOInput

export PWInput,
    namelists,
    cards

@with_kw struct PWInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomic_species::AtomicSpeciesCard
    atomic_positions::AtomicPositionsCard
    k_points::KPointsCard
    cell_parameters::CellParametersCard
end  # struct PWInput

filter_field_by_supertype(obj, ::Type{T}) where {T} = filter(x->isa(x, T), map(x->getfield(obj, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWInput) = filter_field_by_supertype(input, Card)

end
