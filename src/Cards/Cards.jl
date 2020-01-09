"""
# module Cards



# Examples

```jldoctest
julia>
```
"""
module Cards

using QuantumESPRESSOBase: InputEntry

export Card, optionof, allowed_options

abstract type Card <: InputEntry end

"""
    optionof(x::Card)

Return the option for `Card` `x`.

A user should not use `x.option` to access a `Card`'s `option`. Because some `Card`s do not have an option.
Using `optionof(x)` is suggested.
"""
optionof(card::Card) = getfield(card, :option)

"""
    allowed_options(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards, QuantumESPRESSOBase.Cards.PWscf

julia> allowed_options(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> allowed_options(CellParametersCard)
("alat", "bohr", "angstrom")

julia> allowed_options(KPointsCard)
("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
allowed_options(::Type{<:Card}) = nothing

include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")

end
