"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using AbInitioSoftwareBase.Inputs: InputEntry, Namelist, Card, Setter
using AutoHashEquals: @auto_hash_equals
using Compat: isnothing
using OptionalArgChecks: @argcheck
using Setfield: @set!
using Unitful: AbstractQuantity, Length, Temperature, ustrip, @u_str
using UnitfulAtomic

using ..Inputs: QuantumESPRESSOInput, VerbositySetter, entryname

import AbInitioSoftwareBase.Inputs: asstring, groupname
import Crystallography: Bravais, Lattice, cellvolume
# import Pseudopotentials: pseudoformat
import ..Inputs: optionof

export optionof, groupname

include("namelists.jl")
include("cards.jl")
include("input.jl")
include("crystallography.jl")
include("asstring.jl")
include("set.jl")
include("validation.jl")

end
