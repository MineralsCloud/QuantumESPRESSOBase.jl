"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using Setfield: @set!

using QuantumESPRESSOBase: bravais_lattice
using ..Namelists: Namelist
using ..Cards
using ..Cards: Card

export namelists, cards, autofill_cell_parameters, compulsory_namelists, compulsory_cards

abstract type QuantumESPRESSOInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

"""
    namelists(input::QuantumESPRESSOInput)

Return a vector of `Namelist`s of a `QuantumESPRESSOInput`'s subtypes.
"""
namelists(input::QuantumESPRESSOInput) = _filterfields(input, Namelist)

"""
    cards(input::QuantumESPRESSOInput)

Return a vector of `Card`s of a `QuantumESPRESSOInput`'s subtypes.
"""
cards(input::QuantumESPRESSOInput) = _filterfields(input, Card)

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

include("postlude.jl")

end
