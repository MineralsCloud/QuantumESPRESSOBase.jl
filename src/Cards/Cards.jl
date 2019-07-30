"""
# module Cards



# Examples

```jldoctest
julia>
```
"""
module Cards

using Parameters

using QuantumESPRESSOBase: InputEntry

export Card, option, allowed_options

abstract type Card <: InputEntry end

function Parameters.reconstruct(card::Card, newdict::AbstractDict)
    :option in keys && error("If you want to change the option of a card, reconstruct a new one!")
    return reconstruct(card, newdict)
end # function reconstruct

include("PW.jl")
include("option.jl")
include("allowed_options.jl")

end
