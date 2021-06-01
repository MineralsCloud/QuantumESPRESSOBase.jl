"""
    asstring(data::AtomicSpecies)

Return a `String` representing a `AtomicSpecies`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::AtomicSpecies)
    return join(
        (
            indent(data),
            sprintf1("%3s", data.atom),
            sprintf1(floatfmt(data), data.mass),
            data.pseudopot,
        ),
        delimiter(data),
    )
end
"""
    asstring(card::AtomicSpeciesCard)

Return a `String` representing a `AtomicSpeciesCard`, valid for Quantum ESPRESSO's input.
"""
asstring(card::AtomicSpeciesCard) =
    join(("ATOMIC_SPECIES", map(asstring, unique(card.data))...), newline(card))
"""
    asstring(data::AtomicPosition)

Return a `String` representing a `AtomicPosition`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::AtomicPosition)
    content = join(
        (
            indent(data),
            sprintf1("%3s", data.atom),
            map(x -> sprintf1(floatfmt(data), x), data.pos)...,
        ),
        delimiter(data),
    )
    if !all(data.if_pos)
        return join((content, map(Int, data.if_pos)...), delimiter(data))
    else
        return content
    end
end
"""
    asstring(card::AtomicPositionsCard)

Return a `String` representing a `AtomicPositionsCard`, valid for Quantum ESPRESSO's input.
"""
asstring(card::AtomicPositionsCard) = join(
    ("ATOMIC_POSITIONS { $(optionof(card)) }", map(asstring, card.data)...),
    newline(card),
)
"""
    asstring(card::CellParametersCard)

Return a `String` representing a `CellParametersCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::CellParametersCard)
    return join(
        (
            "CELL_PARAMETERS { $(optionof(card)) }",
            map(eachrow(card.data)) do row
                join((sprintf1(floatfmt(card), x) for x in row))
            end...,
        ),
        newline(card),
    )
end
"""
    asstring(data::MonkhorstPackGrid)

Return a `String` representing a `MonkhorstPackGrid`, valid for Quantum ESPRESSO's input.
"""
function asstring(data::MonkhorstPackGrid)
    return indent(data) * join(map([data.mesh; data.is_shift]) do x
        sprintf1(intfmt(data), x)
    end, delimiter(data))
end
"""
    asstring(data::SpecialKPoint)

Return a `String` representing a `SpecialKPoint`, valid for Quantum ESPRESSO's input.
"""
asstring(data::SpecialPoint) =
    indent(data) * join(map(x -> sprintf1(floatfmt(data), x), data), delimiter(data))
"""
    asstring(card::KPointsCard)

Return a `String` representing a `KPointsCard`, valid for Quantum ESPRESSO's input.
"""
function asstring(card::SpecialPointsCard)
    content = "K_POINTS { $(optionof(card)) }" * newline(card)
    return join((content, length(card.data), map(asstring, card.data)...), newline(card))
end
asstring(card::GammaPointCard) = "K_POINTS { $(optionof(card)) }" * newline(card)
function asstring(card::KMeshCard)
    content = "K_POINTS { $(optionof(card)) }" * newline(card)
    return content * asstring(card.data)
end

indent(::Union{AtomicSpecies,AtomicPosition,SpecialPoint,MonkhorstPackGrid,AtomicForce}) =
    ' '^4

delimiter(
    ::Union{AtomicSpecies,AtomicPosition,SpecialPoint,MonkhorstPackGrid,AtomicForce},
) = ' '

floatfmt(::Union{AtomicSpecies,AtomicPosition,SpecialPoint}) = "%14.9f"
floatfmt(::CellParametersCard) = "%14.9f"

intfmt(::MonkhorstPackGrid) = "%5d"
