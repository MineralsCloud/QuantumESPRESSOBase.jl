using Fortran90Namelists.JuliaToFortran: to_fortran
using IterTools: fieldvalues

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs
using QuantumESPRESSOBase.Inputs.PWscf

export to_qe

"""
    to_qe()



# Arguments

# Examples

```jldoctest
julia>
```
"""
function to_qe(dict::AbstractDict; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    content = ""
    f = string ∘ to_fortran
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                ismissing(x) && continue
                content *= "$indent$key($i)$sep=$sep$(f(x))\n"
            end
        else
            content *= "$indent$key$sep=$sep$(f(value))\n"
        end
    end
    return content
end
function to_qe(nml::Namelist; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    namelist_name = (uppercase ∘ string ∘ name ∘ typeof)(nml)
    content = """&$namelist_name
    $(to_qe(dropdefault(nml); indent = indent, sep = sep))/
    """
end
function to_qe(data::AtomicSpecies; sep::AbstractString = " ")::String
    return join(map(string, fieldvalues(data)), sep)
end
function to_qe(card::AtomicSpeciesCard; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    """
    ATOMIC_SPECIES
    $(join(["$(indent)$(to_qe(x; sep = sep))" for x in card.data], "\n"))
    """
end
function to_qe(data::AtomicPosition; sep::AbstractString = " ", with_if_pos::Bool = false)::String
    with_if_pos && return join(map(string, [data.atom; data.pos; data.if_pos]), sep)
    return join(map(string, [data.atom; data.pos]), sep)
end
function to_qe(card::AtomicPositionsCard; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    """
    ATOMIC_POSITIONS$(sep){ $(card.option) }
    $(join(["$(indent)$(to_qe(x; sep = sep))" for x in card.data], "\n"))
    """
end
function to_qe(card::CellParametersCard; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    """
    CELL_PARAMETERS$(sep){ $(card.option) }
    $(join(["$(indent)$(join(row, sep))" for row in eachrow(card.data)], "\n"))
    """
end
function to_qe(data::MonkhorstPackGrid; sep::AbstractString = " ")::String
    return join(map(string, [data.grid; data.offsets]), sep)
end
function to_qe(data::GammaPoint)::String
    return ""
end
function to_qe(data::SpecialKPoint; sep::AbstractString = " ")::String
    return join(map(string, [data.coordinates; data.weight]), sep)
end
function to_qe(card::KPointsCard; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    content = "K_POINTS$(sep){ $(card.option) }\n"
    if card.option in ("gamma", "automatic")
        content *= "$(indent)$(to_qe(first(card.data)))\n"
    else  # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= "$(length(card.data))\n"
        for x in card.data
            content *= "$(indent)$(to_qe(x, sep = sep))\n"
        end
    end
    return content
end
function to_qe(input::PWscfInput; indent::AbstractString = "    ", sep::AbstractString = " ", verbose::Bool = true)::String
    if verbose
        return join(map(to_qe, fieldvalues(input)), "\n")
    else
        str = ""
        for namelist in namelists(input)
            str *= to_qe(namelist)
        end
        for card in cards(input)
            str *= to_qe(card)
        end
        return str
    end
end
