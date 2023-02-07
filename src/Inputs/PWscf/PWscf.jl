"""
# module PWscf



# Examples

```julia
julia>
```
"""
module PWscf

using AbInitioSoftwareBase.Inputs: InputEntry, Namelist, Card, Setter
using Setfield: @set!
using StructEquality: @struct_hash_equal
using Unitful: AbstractQuantity, Length, Temperature, ustrip, @u_str
using UnitfulAtomic

using ..Inputs: QuantumESPRESSOInput, VerbositySetter

import AbInitioSoftwareBase.Inputs: groupname, getpseudodir, getpotentials
import CrystallographyBase: Bravais, Lattice, cellvolume
# import Pseudopotentials: pseudoformat
import ..Inputs: optionof

export optionof, groupname

include("namelists.jl")
include("cards.jl")
include("input.jl")
include("crystallography.jl")
include("set.jl")
include("validation.jl")

end
