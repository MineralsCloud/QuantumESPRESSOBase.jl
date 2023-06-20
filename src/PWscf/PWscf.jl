module PWscf

using AbInitioSoftwareBase: InputEntry, Namelist, Card, Setter
using Setfield: @set!
using StructEquality: @struct_hash_equal
using Unitful: AbstractQuantity, Length, Temperature, ustrip, @u_str
using UnitfulAtomic

using ..QuantumESPRESSOBase: QuantumESPRESSOInput, VerbositySetter

import AbInitioSoftwareBase: groupname, getpseudodir
import CrystallographyBase: Lattice, cellvolume
# import Pseudopotentials: pseudoformat

export groupname

include("namelists.jl")
include("cards.jl")
include("inputs.jl")
include("crystallography.jl")
include("set.jl")
include("validation.jl")

end
