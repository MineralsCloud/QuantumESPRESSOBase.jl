using LinearAlgebra: det

using ..QuantumESPRESSOBase: Ibrav, latticevectors

import CrystallographyBase: Cell, crystaldensity

struct InsufficientInfoError <: Exception
    msg::AbstractString
end

"""
    Lattice(nml::SystemNamelist)

Create a `Lattice` from a `SystemNamelist`.
"""
Lattice(nml::SystemNamelist) = Lattice(latticevectors(nml.celldm, Ibrav(nml.ibrav)))
"""
    Lattice(card::CellParametersCard)

Create a `Lattice` from a `CellParametersCard`.
"""
function Lattice(card::CellParametersCard)
    m, option = transpose(card.data), getoption(card)
    if option == :alat
        throw(InsufficientInfoError("parameter `celldm[1]` needed!"))
    elseif option == :bohr
        return Lattice(m)
    else  # option == :angstrom
        return Lattice(m * ustrip(u"bohr", 1u"angstrom"))
    end
end
"""
    Lattice(card::PWInput)

Create a `Lattice` from a `PWInput`.
"""
function Lattice(input::PWInput)
    if isnothing(input.cell_parameters)
        return Lattice(input.system)
    else
        if getoption(input.cell_parameters) == "alat"
            return Lattice(
                transpose(input.cell_parameters.data) * first(input.system.celldm)
            )
        else
            return Lattice(input.cell_parameters)
        end
    end
end

function Cell(input::PWInput)
    lattice = Lattice(input) * 1u"bohr"
    positions = [position for (_, position) in eachatom(input.atomic_positions)]
    atoms = [atom for (atom, _) in eachatom(input.atomic_positions)]
    return Cell(lattice, positions, atoms)
end

"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function cellvolume(card::AbstractCellParametersCard)
    option = getoption(card)
    if option == :bohr
        return abs(det(card.data))
    elseif option == :angstrom
        return ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == :alat
        throw(InsufficientInfoError("parameter `celldm[1]` needed!"))
    end
end
"""
    cellvolume(nml::SystemNamelist)

Return the volume of the cell based on the information given in a `SystemNamelist`, in atomic unit.
"""
cellvolume(nml::SystemNamelist) = cellvolume(Lattice(nml))
"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function cellvolume(input::PWInput)
    if input.system.ibrav == 0
        if isnothing(input.cell_parameters)
            throw(InsufficientInfoError("`ibrav` is 0, must read cell parameters!"))
        else
            if getoption(input.cell_parameters) == :alat
                # If no value of `celldm` is changed...
                if isnothing(input.system.celldm[1])
                    throw(InsufficientInfoError("parameter `celldm[1]` needed!"))
                else
                    return input.system.celldm[1]^3 * abs(det(input.cell_parameters.data))
                end
            else  # "bohr" or "angstrom"
                return cellvolume(input.cell_parameters)
            end
        end
    else
        return cellvolume(Lattice(input.system))
    end
end

function crystaldensity(input::PWInput)
    lattice = Lattice(input) * 1u"bohr"
    atoms = (Symbol(uppercasefirst(atom)) for (atom, _) in eachatom(input.atomic_positions))
    return crystaldensity(lattice, atoms)
end
