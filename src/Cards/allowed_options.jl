using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf

export allowed_options

"""
    allowed_options(T::Type{<: Card})

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
allowed_options(::Type{<:AtomicPositionsCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{<:CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{<:KPointsCard}) =
    ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
