using CrystallographyBase: CartesianFromFractional
using LinearAlgebra: det
using Spglib: get_dataset

using ..Inputs: Ibrav

import CrystallographyBase: Cell, crystaldensity

export find_symmetry

struct InsufficientInfoError <: Exception
    msg::AbstractString
end

"""
    Ibrav(nml::SystemNamelist)

Return a `Ibrav` from a `SystemNamelist`.
"""
Ibrav(nml::SystemNamelist) = Ibrav(nml.ibrav)

"""
    Lattice(nml::SystemNamelist)

Create a `Lattice` from a `SystemNamelist`.
"""
Lattice(nml::SystemNamelist) = Lattice(nml.celldm, Ibrav(nml))
"""
    Lattice(card::CellParametersCard)

Create a `Lattice` from a `CellParametersCard`.
"""
function Lattice(card::CellParametersCard)
    m, option = transpose(card.data), optionof(card)
    if option == "alat"
        throw(InsufficientInfoError("parameter `celldm[1]` needed!"))
    elseif option == "bohr"
        return Lattice(m)
    else  # option == "angstrom"
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
        if optionof(input.cell_parameters) == "alat"
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
    positions = [atomic_position.pos for atomic_position in input.atomic_positions.data]
    types = [atomic_position.atom for atomic_position in input.atomic_positions.data]
    return Cell(lattice, positions, types)
end

"""
    cellvolume(card)

Return the cell volume of a `CellParametersCard` or `RefCellParametersCard`, in atomic unit.

!!! warning
    It will throw an error if the option is `"alat"`.
"""
function cellvolume(card::AbstractCellParametersCard)
    option = optionof(card)
    if option == "bohr"
        return abs(det(card.data))
    elseif option == "angstrom"
        return ustrip(u"bohr^3", abs(det(card.data)) * u"angstrom^3")
    else  # option == "alat"
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
            if optionof(input.cell_parameters) == "alat"
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
    atoms = (
        Symbol(uppercasefirst(atomic_position.atom)) for
        atomic_position in input.atomic_positions.data
    )
    return crystaldensity(lattice, atoms)
end

function find_symmetry(input::PWInput, symprec=1e-5)
    lattice = Lattice(input)
    option = input.atomic_positions.option
    data = Iterators.map(input.atomic_positions.data) do atomic_position
        atom, position = atomic_position.atom, atomic_position.pos
        # `position` is a `Vector` in unit of "bohr"
        if option == "alat"
            position *= input.system.celldm[1]
        elseif option == "bohr"
            position
        elseif option == "angstrom"
            ustrip.(u"bohr", position * u"angstrom")
        elseif option == "crystal"
            CartesianFromFractional(lattice)(position)
        else  # option == "crystal_sg"
            error("unimplemented!")  # FIXME
        end
        position, atom
    end
    cell = Cell(lattice, first.(data), last.(data))
    dataset = get_dataset(cell, symprec)
    return dataset
end
