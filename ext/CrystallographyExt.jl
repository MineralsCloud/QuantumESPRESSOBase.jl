module CrystallographyExt

using CrystallographyBase: Lattice, Cell
using QuantumESPRESSOBase.PWscf: PWInput
using Unitful: ustrip, @u_str
using UnitfulAtomic

import Crystallography: findsymmetry

function findsymmetry(input::PWInput, symprec=1e-5)
    lattice = Lattice(input)
    option = input.atomic_positions.option
    data = Iterators.map(input.atomic_positions.data) do atomic_position
        atom, position = atomic_position.atom, atomic_position.pos
        # `position` is a `Vector` in unit of "bohr"
        if option == :alat
            position *= input.system.celldm[1]
        elseif option == :bohr
            position
        elseif option == :angstrom
            ustrip.(u"bohr", position .* u"angstrom")
        elseif option == :crystal
            lattice(position)
        else  # option == :crystal_sg
            error("unimplemented!")  # FIXME
        end
        position, atom
    end
    cell = Cell(lattice, first.(data), last.(data))
    return findsymmetry(cell, symprec)
end

end
