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

include("prelude.jl")
include("PWscf.jl")
include("option.jl")
include("allowed_options.jl")

end
