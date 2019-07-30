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
name(::Type{<: InputEntry}) = error("Undefined name!")
name(::Type{<: ControlNamelist}) = :control
name(::Type{<: SystemNamelist}) = :system
name(::Type{<: ElectronsNamelist}) = :electrons
name(::Type{<: IonsNamelist}) = :ions
name(::Type{<: CellNamelist}) = :cell
name(::Type{<: INPUTPHNamelist}) = :inputph
name(::Type{<: AtomicSpeciesCard}) = :atomicspecies
name(::Type{<: AtomicPositionsCard}) = :atomicpositions
name(::Type{<: KPointsCard}) = :kpoints
name(::Type{<: CellParametersCard}) = :cellparameters
