using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf

export option

"""
    option(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `option(x)` is suggested.
"""
option(card::Card) = getfield(card, :option)
option(card::AtomicSpeciesCard) = nothing
