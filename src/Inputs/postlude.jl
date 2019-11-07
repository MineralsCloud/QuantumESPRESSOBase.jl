using .PWscf
using .CP

"""
    autofill_cell_parameters(template::Union{PWscfInput,CPInput})

Generate automatically a `CellParametersCard` for a `PWscfInput` or `CPInput` if its `cell_parameters` field is `nothing`.

Sometimes the `ibrav` field of a `PWscfInput` is not `0`, with its `cell_parameters` field to be empty.
But there are cases we want to write its `CellParametersCard` explicitly. This function will take a `PWscfInput` described
above and generate a new `PWscfInput` with its `ibrav = 0` and `cell_parameters` not empty.
"""
function autofill_cell_parameters(template::Union{PWscfInput,CPInput})
    system = template.system
    @set! template.cell_parameters = CellParametersCard("alat", bravais_lattice(system))
    @set! template.system.ibrav = 0
    @set! template.system.celldm = [system.celldm[1]]
end # function autofill_cell_parameters

"""
    compulsory_namelists(input::Union{PWscfInput,CPInput})

Return a vector of compulsory `Namelist`s of a `PWscfInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
"""
compulsory_namelists(input::Union{PWscfInput,CPInput}) = [getfield(input, x) for x in (:control, :system, :electrons)]

"""
    compulsory_cards(input::PWscfInput)

Return a vector of compulsory `Card`s of a `PWscfInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
"""
compulsory_cards(input::PWscfInput) = [getfield(input, x) for x in (:atomic_species, :atomic_positions, :k_points)]
"""
    compulsory_cards(input::CPInput)

Return a vector of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
"""
compulsory_cards(input::CPInput) = [getfield(input, x) for x in (:atomic_species, :atomic_positions)]

function Cards.cell_volume(input::PWscfInput)
    if isnothing(input.cell_parameters)
        return det(bravais_lattice(input.system))
    else
        iszero(input.system.celldm[1]) && return cell_volume(input.cell_parameters)
        return input.system.celldm[1]^3 * cell_volume(input.cell_parameters)
    end
end # function cell_volume
