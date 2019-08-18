using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf

export option

"""
    option()



# Arguments

# Examples

```jldoctest
julia>
```
"""
option(card::Card) = getfield(card, :option)
option(card::AtomicSpeciesCard) = nothing
