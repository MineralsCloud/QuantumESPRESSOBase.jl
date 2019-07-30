using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Namelists.PH
using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput
using QuantumESPRESSO.QuantumESPRESSOInput.PW

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
