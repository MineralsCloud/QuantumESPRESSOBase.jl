using LinearAlgebra: det
using Unitful: uparse, unitmodules

using ..Inputs: Ibrav

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
    m, option = transpose(card.data), card.option
    return if option == "alat"
        Lattice(m)
    else  # option in ("bohr", "angstrom")
        Lattice(m * uparse(option; unit_context = unitmodules))
    end
end

struct InformationNotEnough <: Exception
    msg::AbstractString
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
        throw(InformationNotEnough("parameter `celldm[1]` needed!"))
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
            throw(InformationNotEnough("`ibrav` is 0, must read cell parameters!"))
        else
            if optionof(input.cell_parameters) == "alat"
                # If no value of `celldm` is changed...
                if isnothing(input.system.celldm[1])
                    throw(InformationNotEnough("parameter `celldm[1]` needed!"))
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
