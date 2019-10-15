using .Namelists
using .Namelists.PWscf
using .Namelists.PHonon
using .Cards
using .Cards.PWscf
using .Inputs
using .Inputs.PWscf

export asfieldname

"""
    asfieldname()



# Arguments

# Examples

```jldoctest
julia>
```
"""
asfieldname(::Type{T}) where {T<:InputEntry} = error("Undefined for entry $(nameof(T))!")
asfieldname(::Type{<:ControlNamelist}) = :control
asfieldname(::Type{<:SystemNamelist}) = :system
asfieldname(::Type{<:ElectronsNamelist}) = :electrons
asfieldname(::Type{<:IonsNamelist}) = :ions
asfieldname(::Type{<:CellNamelist}) = :cell
asfieldname(::Type{<:PHNamelist}) = :inputph
asfieldname(::Type{<:Q2RNamelist}) = :input
asfieldname(::Type{<:MatdynNamelist}) = :input
asfieldname(::Type{<:DynmatNamelist}) = :input
asfieldname(::Type{<:AtomicSpeciesCard}) = :atomic_species
asfieldname(::Type{<:AtomicPositionsCard}) = :atomic_positions
asfieldname(::Type{<:KPointsCard}) = :k_points
asfieldname(::Type{<:CellParametersCard}) = :cell_parameters
