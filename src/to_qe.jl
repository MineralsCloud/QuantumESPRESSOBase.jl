using Fortran90Namelists.JuliaToFortran: to_fortran
using IterTools: fieldvalues

using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput
using QuantumESPRESSO.QuantumESPRESSOInput.PW

export to_qe

"""
    to_qe()



# Arguments

# Examples

```jldoctest
julia>
```
"""
function to_qe(dict::AbstractDict; indent::AbstractString = "    ")::String
    content = ""
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                ismissing(x) && continue
                content *= "$(indent)$(key)($i) = $(string(to_fortran(x)))\n"
            end
        else
            content *= "$(indent)$(key) = $(string(to_fortran(value)))\n"
        end
    end
    return content
end
function to_qe(nml::Namelist; indent::AbstractString = "    ")::String
    return "&$(name(typeof(nml)))\n" * to_qe(to_dict(nml); indent = indent) * "/\n"
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
function to_qe(input::PWInput; indent::AbstractString = "    ", sep::AbstractString = " ", debug::Bool = true)::String
    if debug
        return join(map(to_qe, fieldvalues(input)), "\n")
    else
        str = ""
        for namelist in namelists(input)
            str *= "&" * uppercase(string(name(typeof(namelist)))) * "\n"
            str *= to_qe(dropdefault(namelist))
            str *= "/\n"
        end
        for card in cards(input)
            str *= to_qe(card)
        end
        return str
    end
end
