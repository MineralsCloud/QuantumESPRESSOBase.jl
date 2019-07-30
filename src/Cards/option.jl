using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW

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
