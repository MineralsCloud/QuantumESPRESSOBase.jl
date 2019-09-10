#=
NamelistsTests.jl:
- Julia version: 1.0
- Author: qz
- Date: Jul 11, 2019
=#
module NamelistsTests

using Test

using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs.PWscf

as = AtomicSpeciesCard([AtomicSpecies("Fe", 55.845, "Fe.pseudopotential")])
ap = AtomicPositionsCard(data = [AtomicPosition(atom = "Fe", pos = [0, 0, 0])])
cell = CellParametersCard(data = ones(3, 3))
k = KPointsCard(option = "gamma", data = [GammaPoint()])

pw = PWscfInput(atomic_species = as, atomic_positions = ap, k_points = k, cell_parameters = cell)

end
