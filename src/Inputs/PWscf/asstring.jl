using AbInitioSoftwareBase.Inputs: FormatConfig
using Compat: eachrow
using Formatting: sprintf1

export asstring

formatconfig(
    ::Type{
        <:Union{
            AtomicSpecies,
            AtomicPosition,
            ReciprocalPoint,
            AtomicSpeciesCard,
            CellParametersCard,
            MonkhorstPackGrid,
            AtomicForce,
        },
    },
) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%14.9f",
    int = "%5d",
    bool = ".%.",
)

"""
    asstring(data::AtomicSpecies)

Return a `String` representing a `AtomicSpecies`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::AtomicSpecies)
    config = formatconfig(data)
    return join(
        (
            config.indent,
            sprintf1("%3s", data.atom),
            sprintf1(config.float, data.mass),
            data.pseudopot,
        ),
        config.delimiter,
    )
end
"""
    asstring(card::AtomicSpeciesCard)

Return a `String` representing a `AtomicSpeciesCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::AtomicSpeciesCard)
    config = formatconfig(card)
    data = union(card.data)
    return join(("ATOMIC_SPECIES", map(asstring, data)...), config.newline)
end
"""
    asstring(data::AtomicPosition)

Return a `String` representing a `AtomicPosition`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::AtomicPosition)
    config = formatconfig(data)
    content = join(
        (
            config.indent,
            sprintf1("%3s", data.atom),
            map(x -> sprintf1(config.float, x), data.pos)...,
        ),
        config.delimiter,
    )
    if !all(data.if_pos)
        return join((content, map(Int, data.if_pos)...), config.delimiter)
    else
        return content
    end
end
"""
    asstring(card::AtomicPositionsCard)

Return a `String` representing a `AtomicPositionsCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::AtomicPositionsCard)
    config = formatconfig(card)
    join(
        ("ATOMIC_POSITIONS { $(optionof(card)) }", map(asstring, card.data)...),
        config.newline,
    )
end
"""
    asstring(card::CellParametersCard)

Return a `String` representing a `CellParametersCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::CellParametersCard)
    config = formatconfig(card)
    return join(
        (
            "CELL_PARAMETERS { $(optionof(card)) }",
            map(eachrow(card.data)) do row
                join((sprintf1(config.float, x) for x in row))
            end...,
        ),
        config.newline,
    )
end
"""
    asstring(data::MonkhorstPackGrid)

Return a `String` representing a `MonkhorstPackGrid`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::MonkhorstPackGrid)
    config = formatconfig(data)
    return config.indent * join(map([data.mesh; data.is_shift]) do x
        sprintf1(config.int, x)
    end, config.delimiter)
end
"""
    asstring(data::SpecialKPoint)

Return a `String` representing a `SpecialKPoint`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::ReciprocalPoint)
    config = formatconfig(data)
    return config.indent * join(
        map(x -> sprintf1(config.float, x), [data.coord..., data.weight]),
        config.delimiter,
    )
end
"""
    asstring(card::KPointsCard)

Return a `String` representing a `KPointsCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::SpecialPointsCard)
    config = formatconfig(card)
    content = "K_POINTS { $(optionof(card)) }" * config.newline
    return join((content, length(card.data), map(asstring, card.data)...), config.newline)
end
function asstring(card::GammaPointCard)
    config = formatconfig(card)
    return "K_POINTS { $(optionof(card)) }" * config.newline
end
function asstring(card::KMeshCard)
    config = formatconfig(card)
    content = "K_POINTS { $(optionof(card)) }" * config.newline
    return content * asstring(card.data)
end
