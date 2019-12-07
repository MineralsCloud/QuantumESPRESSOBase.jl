"""
# module Inputs



# Examples

```jldoctest
julia>
```
"""
module Inputs

using LinearAlgebra: det

using Compat: isnothing
using Setfield: @set!

using QuantumESPRESSOBase
using QuantumESPRESSOBase: bravais_lattice
using ..Namelists: Namelist
using ..Cards
using ..Cards: Card, optionof

export namelists, cards, autofill_cell_parameters, compulsory_namelists, compulsory_cards

abstract type QuantumESPRESSOInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(f, obj) = Iterators.filter(f, (getfield(obj, i) for i in 1:nfields(obj)))

"""
    namelists(input::QuantumESPRESSOInput)

Return an iterable of `Namelist`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
namelists(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Namelist), input)

"""
    cards(input::QuantumESPRESSOInput)

Return an iterable of `Card`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
cards(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Card), input)

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

using .PWscf: PWInput
using .CP: CPInput

"""
    autofill_cell_parameters(template::Union{PWInput,CPInput})

Generate automatically a `CellParametersCard` for a `PWInput` or `CPInput` if its `cell_parameters` field is `nothing`.

Sometimes the `ibrav` field of a `PWInput` is not `0`, with its `cell_parameters` field to be empty.
But there are cases we want to write its `CellParametersCard` explicitly. This function will take a `PWInput` described
above and generate a new `PWInput` with its `ibrav = 0` and `cell_parameters` not empty.
"""
function autofill_cell_parameters(template::Union{PWInput,CPInput})
    system = template.system
    @set! template.cell_parameters =
        Cards.CellParametersCard("alat", bravais_lattice(system))
    @set! template.system.ibrav = 0
    @set! template.system.celldm = [system.celldm[1]]
end # function autofill_cell_parameters

"""
    compulsory_namelists(input::Union{PWInput,CPInput})

Return an iterable of compulsory `Namelist`s of a `PWInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
It is lazy, you may want to `collect` it.

To get the optional `Namelist`s, use `(!compulsory_namelists)(input)`.
"""
compulsory_namelists(input::Union{PWInput,CPInput}) =
    (getfield(input, x) for x in (:control, :system, :electrons))
Base.:!(::typeof(compulsory_namelists)) = input -> setdiff(collect(namelists(input)), collect(compulsory_namelists(input)))

"""
    compulsory_cards(input::PWInput)

Return an iterable of compulsory `Card`s of a `PWInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)`.
"""
compulsory_cards(input::PWInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions, :k_points))
"""
    compulsory_cards(input::CPInput)

Return an iterable of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)`.
"""
compulsory_cards(input::CPInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions))
Base.:!(::typeof(compulsory_cards)) = input -> setdiff(collect(cards(input)), collect(compulsory_cards(input)))

function QuantumESPRESSOBase.cell_volume(input::PWInput)
    if isnothing(input.cell_parameters)
        return det(bravais_lattice(input.system))
    else
        if optionof(input.cell_parameters) == "alat"
            # If no value of `celldm` is changed...
            isnothing(input.system.celldm[1]) && error("`celldm[1]` is not defined!")
            return input.system.celldm[1]^3 * det(input.cell_parameters.data)
        else  # "bohr" or "angstrom"
            return cell_volume(input.cell_parameters)
        end
    end
end # function QuantumESPRESSOBase.cell_volume

function QuantumESPRESSOBase.to_qe(
    input::QuantumESPRESSOInput;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
    verbose::Bool = false,
)::String
    content = ""
    for namelist in namelists(input)
        content *= to_qe(namelist, indent = indent, sep = sep, verbose = verbose)
    end
    for card in cards(input)
        content *= to_qe(card, indent = indent, sep = sep)
    end
    return content
end

end
