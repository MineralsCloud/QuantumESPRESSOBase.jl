using Crystallography: CartesianFromFractional
using LinearAlgebra: det
using Spglib: get_dataset

using ..Inputs: Ibrav

export find_symmetry

struct LackCellInfoError <: Exception
    msg::AbstractString
end

"""
    Bravais(nml::SystemNamelist)

Return a `Bravais` from a `SystemNamelist`.
"""
Bravais(nml::SystemNamelist) = Bravais(Ibrav(nml.ibrav))

"""
    Lattice(nml::SystemNamelist)

Return a `Lattice` from a `SystemNamelist`.
"""
Lattice(nml::SystemNamelist) = Lattice(Bravais(nml), nml.celldm)
"""
    Lattice(card::CellParametersCard)

Return a `Lattice` from a `CellParametersCard`.
"""
function Lattice(card::CellParametersCard)
    m, option = transpose(card.data), optionof(card)
    if option == "alat"
        throw(LackCellInfoError("parameter `celldm[1]` needed!"))
    elseif option == "bohr"
        return Lattice(m)
    else  # option == "angstrom"
        return Lattice(m * ustrip(u"bohr", 1u"angstrom"))
    end
end
"""
    Lattice(card::PWInput)

Return a `Lattice` from a `PWInput`.
"""
function Lattice(input::PWInput)
    if isnothing(input.cell_parameters)
        return Lattice(input.system)
    else
        if optionof(input.cell_parameters) == "alat"
            return Lattice(transpose(input.cell_parameters.data) * first(input.system.celldm))
        else
            return Lattice(input.cell_parameters)
        end
    end
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
        throw(LackCellInfoError("parameter `celldm[1]` needed!"))
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
            throw(LackCellInfoError("`ibrav` is 0, must read cell parameters!"))
        else
            if optionof(input.cell_parameters) == "alat"
                # If no value of `celldm` is changed...
                if isnothing(input.system.celldm[1])
                    throw(LackCellInfoError("parameter `celldm[1]` needed!"))
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

function find_symmetry(input::PWInput, symprec = 1e-5)
    lattice = Lattice(input)
    option = input.atomic_positions.option
    data = map(input.atomic_positions.data) do atomic_position
        atom, position = atomic_position.atom, atomic_position.pos
        # In unit of bohr
        if option == "alat"
            position *= input.system.celldm[1]
        elseif option == "bohr"
            position
        elseif option == "angstrom"
            ustrip(u"bohr", position * u"angstrom")
        elseif option == "crystal"
            CartesianFromFractional(lattice)(position)
        else  # option == "crystal_sg"
            error("unimplemented!")  # FIXME
        end
        atom, position
    end
    cell = Cell(lattice.data, first.(data), last.(data))
    dataset = get_dataset(cell, symprec)
    return dataset
end
