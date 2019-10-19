"""
# module Cards



# Examples

```jldoctest
julia>
```
"""
module Cards

using MLStyle: @match
using Parameters: @with_kw
using Setfield: @lens, get

using QuantumESPRESSOBase: InputEntry

export Card, option, allowed_options

include("prelude.jl")
include("common.jl")
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
include("option.jl")
include("allowed_options.jl")

end
