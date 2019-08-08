using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PW
using QuantumESPRESSOBase.Namelists.PH
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PW
using QuantumESPRESSOBase.QuantumESPRESSOInput
using QuantumESPRESSOBase.QuantumESPRESSOInput.PW

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
name(::Type{<: INPUTPHNamelist}) = :inputph
name(::Type{<: AtomicSpeciesCard}) = :atomic_species
name(::Type{<: AtomicPositionsCard}) = :atomic_positions
name(::Type{<: KPointsCard}) = :k_points
name(::Type{<: CellParametersCard}) = :cell_parameters
