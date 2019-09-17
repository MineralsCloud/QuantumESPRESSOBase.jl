"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

export AbstractInput

abstract type AbstractInput end

include("PWscf.jl")
include("PHonon.jl")

end
