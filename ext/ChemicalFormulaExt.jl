module ChemicalFormulaExt

using ChemicalFormula: sumformula
using QuantumESPRESSOBase.PWscf: AtomicPositionsCard, PWInput

import ChemicalFormula: Formula

function Formula(card::AtomicPositionsCard)
    atoms = map(card.data) do position
        filter(isletter, position.atom)
    end
    str = join(symbol^count(atom == symbol for atom in atoms) for symbol in unique(atoms))
    return Formula(sumformula(Formula(str)))
end
Formula(input::PWInput) = Formula(input.atomic_positions)

end
