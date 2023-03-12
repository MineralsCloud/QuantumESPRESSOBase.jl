module PWscf

using AbInitioSoftwareBase.Inputs: InputEntry, Namelist, Card, Setter
using Setfield: @set!
using StructEquality: @struct_hash_equal
using Unitful: AbstractQuantity, Length, Temperature, ustrip, @u_str
using UnitfulAtomic

using ..QuantumESPRESSOBase: QuantumESPRESSOInput, VerbositySetter

import AbInitioSoftwareBase.Inputs: groupname, getpseudodir, getpotentials
import CrystallographyBase: Lattice, cellvolume
# import Pseudopotentials: pseudoformat

export groupname

include("namelists.jl")
include("cards.jl")
include("input.jl")
include("crystallography.jl")
include("set.jl")
include("validation.jl")

end
