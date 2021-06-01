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
using Unitful: AbstractQuantity, Temperature, ustrip, @u_str
using UnitfulAtomic

using ..Inputs: QuantumESPRESSOInput, VerbositySetter, entryname

import AbInitioSoftwareBase.Inputs: asstring, groupname
import Crystallography: Bravais, Lattice, cellvolume
# import Pseudopotentials: pseudoformat
import ..Inputs: optionof

export optionof, groupname

include("namelists.jl")
include("cards.jl")

function iscompatible(system::SystemNamelist, cell_parameters::CellParametersCard)
    ibrav, celldm = system.ibrav, system.celldm
    if iszero(ibrav)
        if optionof(cell_parameters) in ("bohr", "angstrom")
            return all(iszero, celldm)
        else  # "alat"
            return !iszero(first(celldm))  # first(celldm) != 0
        end
    else
        return false
    end
end # function iscompatible
iscompatible(x::CellParametersCard, y::SystemNamelist) = iscompatible(y, x)

include("input.jl")
include("crystallography.jl")
include("asstring.jl")
include("set.jl")

end
