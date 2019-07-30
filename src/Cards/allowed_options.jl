using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW

export allowed_options

"""
    allowed_options()



# Arguments

# Examples

```jldoctest
julia>
```
"""
allowed_options(::Type{<: Card}) = nothing
allowed_options(::Type{<: AtomicPositionsCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{<: CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{<: KPointsCard}) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
