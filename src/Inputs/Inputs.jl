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
using Kaleido: @batchlens

using QuantumESPRESSOBase: BravaisLattice, to_qe
using QuantumESPRESSOBase.Namelists: Namelist
using QuantumESPRESSOBase.Cards: Card, optionof
using QuantumESPRESSOBase.Cards.PWscf: CellParametersCard
using QuantumESPRESSOBase.Setters: CellParametersSetter, LensMaker

import Crystallography.Crystals
import QuantumESPRESSOBase
import QuantumESPRESSOBase.Setters

export namelistsof, cardsof, compulsory_namelists, compulsory_cards

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(f, obj) = Iterators.filter(f, (getfield(obj, i) for i in 1:nfields(obj)))

"""
    namelistsof(input::QuantumESPRESSOInput)

Return an iterable of `Namelist`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
namelistsof(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Namelist), input)

"""
    cardsof(input::QuantumESPRESSOInput)

Return an iterable of `Card`s of a `QuantumESPRESSOInput`. It is lazy, you may want to `collect` it.
"""
cardsof(input::QuantumESPRESSOInput) = _filterfields(x -> isa(x, Card), input)

# =============================== Modules ============================== #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

using .PWscf: PWInput
using .CP: CPInput

"""
    compulsory_namelists(input::Union{PWInput,CPInput})

Return an iterable of compulsory `Namelist`s of a `PWInput` or `CPInput` (`ControlNamelist`, `SystemNamelist` and `ElectronsNamelist`).
It is lazy, you may want to `collect` it.

To get the optional `Namelist`s, use `(!compulsory_namelists)(input)` (Note the parenthesis!).
"""
compulsory_namelists(input::Union{PWInput,CPInput}) =
    (getfield(input, x) for x in (:control, :system, :electrons))
Base.:!(::typeof(compulsory_namelists)) =
    function (input::T) where {T<:Union{PWInput,CPInput}}
        (
            getfield(input, x) for x in fieldnames(T) if x ∉ (
                :control,
                :system,
                :electrons,
            ) && fieldtype(T, x) <: Namelist
        )
    end

"""
    compulsory_cards(input::PWInput)

Return an iterable of compulsory `Card`s of a `PWInput` (`AtomicSpeciesCard`, `AtomicPositionsCard` and `KPointsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)` (Note the parenthesis!).
"""
compulsory_cards(input::PWInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions, :k_points))
"""
    compulsory_cards(input::CPInput)

Return an iterable of compulsory `Card`s of a `CPInput` (`AtomicSpeciesCard` and `AtomicPositionsCard`).
It is lazy, you may want to `collect` it.

To get the optional `Card`s, use `(!compulsory_cards)(input)` (Note the parenthesis!).
"""
compulsory_cards(input::CPInput) =
    (getfield(input, x) for x in (:atomic_species, :atomic_positions))
Base.:!(::typeof(compulsory_cards)) = function (input::T) where {T<:Union{PWInput,CPInput}}
    (
        getfield(input, x) for x in fieldnames(T) if x ∉ (
            :atomic_species,
            :atomic_positions,
        ) && fieldtype(T, x) <: Card
    )
end

"""
    cellvolume(input::PWInput)

Return the volume of the cell based on the information given in a `PWInput`, in atomic unit.
"""
function Crystals.cellvolume(input::PWInput)
    if isnothing(input.cell_parameters)
        return abs(det(BravaisLattice(input.system)()))
    else
        if optionof(input.cell_parameters) == "alat"
            # If no value of `celldm` is changed...
            isnothing(input.system.celldm[1]) && error("`celldm[1]` is not defined!")
            return input.system.celldm[1]^3 * abs(det(input.cell_parameters.data))
        else  # "bohr" or "angstrom"
            return cellvolume(input.cell_parameters)
        end
    end
end # function Crystals.cellvolume

function QuantumESPRESSOBase.to_qe(
    input::QuantumESPRESSOInput;
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
    verbose::Bool = false,
)::String
    content = ""
    for namelist in namelistsof(input)
        content *=
            to_qe(namelist, indent = indent, delim = delim, newline = newline) * newline
    end
    for card in cardsof(input)
        content *= to_qe(card, indent = indent, delim = delim, newline = newline) * newline
    end
    return content
end

function Setters.make(::LensMaker{CellParametersSetter,<:Union{PWInput,CPInput}})
    return @batchlens begin
        _.cell_parameters
        _.system.ibrav
        _.system.celldm
    end
end # function Setters.make

function Setters.preset_values(::CellParametersSetter, template::Union{PWInput,CPInput})
    # !isnothing(template.cell_parameters) && return template
    system = template.system
    return (
        CellParametersCard("alat", BravaisLattice(system)()),
        0,
        [system.celldm[1]],
    )
end # function Setters.preset_values

end
