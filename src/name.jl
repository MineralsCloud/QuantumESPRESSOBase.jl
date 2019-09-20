using .Namelists
using .Namelists.PWscf
using .Namelists.PHonon
using .Cards
using .Cards.PWscf
using .Inputs
using .Inputs.PWscf

export name

"""
    name()



# Arguments

# Examples

```jldoctest
julia>
```
"""
name(::Type{T}) where {T<:InputEntry} = error("Undefined for entry $(nameof(T))!")
name(::Type{<:ControlNamelist}) = :control
name(::Type{<:SystemNamelist}) = :system
name(::Type{<:ElectronsNamelist}) = :electrons
name(::Type{<:IonsNamelist}) = :ions
name(::Type{<:CellNamelist}) = :cell
name(::Type{<:PHNamelist}) = :inputph
name(::Type{<:Q2RNamelist}) = :input
name(::Type{<:MatdynNamelist}) = :input
name(::Type{<:DynmatNamelist}) = :input
name(::Type{<:AtomicSpeciesCard}) = :atomic_species
name(::Type{<:AtomicPositionsCard}) = :atomic_positions
name(::Type{<:KPointsCard}) = :k_points
name(::Type{<:CellParametersCard}) = :cell_parameters
