using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs
using QuantumESPRESSOBase.Inputs.PWscf

export name

"""
    name()



# Arguments

# Examples

```jldoctest
julia>
```
"""
name(::Type{T}) where {T <: InputEntry} = error("Undefined for entry $(nameof(T))!")
name(::Type{<: ControlNamelist}) = :control
name(::Type{<: SystemNamelist}) = :system
name(::Type{<: ElectronsNamelist}) = :electrons
name(::Type{<: IonsNamelist}) = :ions
name(::Type{<: CellNamelist}) = :cell
name(::Type{<: PhononNamelist}) = :inputph
name(::Type{<: AtomicSpeciesCard}) = :atomic_species
name(::Type{<: AtomicPositionsCard}) = :atomic_positions
name(::Type{<: KPointsCard}) = :k_points
name(::Type{<: CellParametersCard}) = :cell_parameters
